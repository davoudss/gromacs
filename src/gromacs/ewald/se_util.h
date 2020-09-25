#ifndef __SE_UTIL_H__
#define __SE_UTIL_H__

#include <stdio.h>
#include <string.h>
#include <cmath>
#include "gromacs/simd/scalar/scalar_math.h"
#include "pme-internal.h"
#include "emmintrin.h"
#include "immintrin.h"
#include "x86intrin.h"

// --------------------------------------------------------------------------
static inline int half(int p)
{
  return (is_odd(p) ? (p-1)/2 : p/2);
}

// -----------------------------------------------------------------------------
static inline int SE_prod3(const int v[3])
{
  return v[0]*v[1]*v[2];
}


// -----------------------------------------------------------------------------
// Set array elements to real-precision zero
static inline void
SE_fp_set_zero(real* x, const int N)
{
  memset(x,0.0,N*sizeof(real));
}


// -----------------------------------------------------------------------------
// unpacking params
static inline void
parse_params(SE_opt* opt, real xi)
{
  real      h0 = opt->box[0]/ (real) opt->M;
  real      w0 = (real) opt->P*h0/2.;
  real    eta0 = (2.0*w0*xi/opt->m)*(2.0*w0*xi/opt->m);
  real      c0 = 2.0*xi*xi/eta0;

  opt->h = h0;
  opt->c = c0;
  opt->w = w0;
  opt->eta = eta0;
  opt->xi = xi;

}

// -----------------------------------------------------------------------------
// packing SE parameters
static inline void
SE_FCN_params(SE_params* params, const SE_opt* opt, int N)
{
  params->N = N;
  params->P = (int) opt->P;
  params->P_half=half( opt->P );
  params->c = opt->c;
  params->d = pow(params->c/PI,1.5);
  params->h = opt->box[0]/opt->M;
  params->a = -FGG_INF;
  params->eta = opt->eta;
  params->beta = opt->beta;
  params->box[XX] = opt->box[XX];
  params->box[YY] = opt->box[YY];
  params->box[ZZ] = opt->box[ZZ];

  params->dims[0] = opt->M;
  params->dims[1] = opt->M;
  params->dims[2] = opt->M;

  params->npdims[0] = params->dims[0]+params->P;
  params->npdims[1] = params->dims[1]+params->P;
  params->npdims[2] = params->dims[2]+params->P;

}

// -----------------------------------------------------------------------------
static void
SE_FGG_base_gaussian(real* zs, const SE_params* params)
{
  int idx ,i,j,k;
  real ih2, ijh2, ijkh2;

  // unpack parameters
  const int p=params->P;
  const int p_half = half(p);
  const int p_from = (is_odd(p) ? p_half:p_half-1);
  const real h2=(params->h)*(params->h);
  const real c=params->c;
  const real d=params->d;
  
  for(i = -p_from; i<=p_half; i++)
    {
      // hoisting this index calculation (more) breaks omp-parallel code
      idx = __IDX3_RMAJ(i+p_from, 0, 0, p, p);   
      ih2 = i*i*h2;
      for(j = -p_from; j<=p_half; j++)
	{
	  ijh2 = ih2+ j*j*h2;
	  for(k = -p_from; k<=p_half; k++)
	    {
	      ijkh2 = ijh2 + k*k*h2;
	      zs[idx++] = d*gmx::exp(-c*ijkh2);
	    }
	}
    }
}



// -----------------------------------------------------------------------------
static int
fgg_expansion(const real x[3], const real q,
	      const SE_params* params,
	      real z2_0[P_MAX],
	      real z2_1[P_MAX],
	      real z2_2[P_MAX],
	      real zf_0[P_MAX],
	      real zf_1[P_MAX],
	      real zf_2[P_MAX])
{
  // unpack params
  const int  p      = params->P;
  const int  p_half = params->P_half;
  const real h      = params->h;
  const real c      = params->c;

  real t0[3];
  int i,j,idx;
  //int idx_from[3]; idx_from is calculated in calc_interpolation_idx 
  // which compute local idx in for parallel cores 
  int p_from;

  // compute index range and centering
  if(is_odd(p))
    {
      for(j=0; j<3; j++)
        {
          idx = (int) ceil(x[j]/h);
          t0[j] = x[j]-h*(idx);
        }
    }
  else
    {
      for(j=0; j<3; j++)
        {
          idx = (int) floor(x[j]/h);
          t0[j] = x[j]-h*(idx);
        }
    }

  // compute third factor 
  real z3 = std::exp(-c*(t0[0]*t0[0] + t0[1]*t0[1] + t0[2]*t0[2]) )*q;

  // compute second factor by induction
  real z_base0 = std::exp(2*c*h*t0[0]);
  real z_base1 = std::exp(2*c*h*t0[1]);
  real z_base2 = std::exp(2*c*h*t0[2]);

  real z0, z1, z2;
  if(is_odd(p))
    {
      z0 = pow(z_base0,-p_half);
      z1 = pow(z_base1,-p_half);
      z2 = pow(z_base2,-p_half);
      p_from = -p_half;
    }
  else
    {
      z0 = pow(z_base0,-p_half+1);
      z1 = pow(z_base1,-p_half+1);
      z2 = pow(z_base2,-p_half+1);
      p_from = -p_half+1;
    }

  z2_0[0] = z0;
  z2_1[0] = z1;
  z2_2[0] = z2;

  // extra terms multiplied to calculate forces
  zf_0[0] = -c*(t0[0]-p_from*h);
  zf_1[0] = -c*(t0[1]-p_from*h);
  zf_2[0] = -c*(t0[2]-p_from*h);

  for(i=1; i<p; i++)
    {
      z0 *=z_base0;
      z1 *=z_base1;
      z2 *=z_base2;

      z2_0[i] = z0;
      z2_1[i] = z1;
      z2_2[i] = z2;

      zf_0[i] = zf_0[i-1]+c*h;
      zf_1[i] = zf_1[i-1]+c*h;
      zf_2[i] = zf_2[i-1]+c*h;
    }

  // save some flops by multiplying one vector with z3 factor
  for(i=0; i<p; i++)
    {
      z2_0[i] *= z3;
    }

  return 0;
}


// ========================================
static real
kaiser(real x, real ow2, real beta) {
  real v = x*x*ow2;
  real t = sqrt(1. - v);
  real e = std::exp(beta*(t-1));
  e = (v<=1) ? e : 0;
  return e;
}

static real
kaiserExtended(real x, real ow2, real beta) {
  real v = x*x*ow2;
  real t = sqrt(1. - v);
  real e = 2.0*std::cosh(beta*t)/std::exp(beta);
  e = (v<=1) ? e : 0;
  return e;
}


static real
dkaiser(real c, real w, real x )
{
  real denom = sqrt(w*w-x*x);
  real dk = c*x/denom;
  dk = (std::abs(denom)<1e-12) ? 0 : dk;
  return dk;
}

static real
dkaiserExtended(real beta, real w, real x )
{
  real v = x/w;
  real y = sqrt(1-v*v);
  real extra_fac = kaiserExtended(x, 1.0/(w*w), beta);
  real dk = -1*beta*v*std::sinh(beta*y)/(y*extra_fac*std::exp(beta))/w;
  dk = (std::abs(y)<1e-12) ? 0 : dk;
  return dk;
}

static
int kaiser_expansion(const real x[3], const real q,
		     const SE_params* params,
		     real z2_0[P_MAX],
		     real z2_1[P_MAX],
		     real z2_2[P_MAX],
		     real zf_0[P_MAX],
		     real zf_1[P_MAX],
		     real zf_2[P_MAX])
{
    // unpack params
    const int p      = params->P;
    const real h     = params->h;
    const real one_h = 1./h;
    const real w     = params->P/2.;
    const real ow2   = 1./(w*w);
    const real beta  = params->beta;
    real t0[3];

    int idx;
    int idx_from;

    int p_half = w;
    
    // compute index range and centering
    if(is_odd(p)) {
      for(int j=0; j<3; j++)
	{
	  idx = (int) ceil(x[j]*one_h);
	  idx_from = idx - p_half;
	  t0[j] = x[j]*one_h-idx_from;
	}
    }
    else {
      for(int j=0; j<3; j++)
	{
	  idx = (int) floor(x[j]*one_h);
	  idx_from = idx - (p_half-1);
	  t0[j] = x[j]*one_h-idx_from;
	}
    }

    // compute second factor
#if 1
    for(int i=0; i<p; i++) {
      z2_0[i] = kaiserExtended(t0[0]-i,ow2,beta);
      z2_1[i] = kaiserExtended(t0[1]-i,ow2,beta);
      z2_2[i] = kaiserExtended(t0[2]-i,ow2,beta);
    }
#else
    for(int i=0; i<p; i++) {
      z2_0[i] = kaiserExtended(t0[0]-i,ow2,beta);
      z2_1[i] = kaiserExtended(t0[1]-i,ow2,beta);
      z2_2[i] = kaiserExtended(t0[2]-i,ow2,beta);
    }
#endif
    
    // save some flops by multiplying one vector with q
    for(int i=0; i<p; i++)
      z2_0[i] *= q;
#if 1
    real c = -beta / (2.*w);
    for (int i=0; i<p; i++){
      zf_0[i] = dkaiser(c, w, t0[0]-i);
      zf_1[i] = dkaiser(c, w, t0[1]-i);
      zf_2[i] = dkaiser(c, w, t0[2]-i);
    }
#else
    for (int i=0; i<p; i++){
      zf_0[i] = dkaiserExtended(beta, w, t0[0]-i);
      zf_1[i] = dkaiserExtended(beta, w, t0[1]-i);
      zf_2[i] = dkaiserExtended(beta, w, t0[2]-i);
    }
#endif

    
    return 0;
}

#endif //__SE_UTIL_H__
