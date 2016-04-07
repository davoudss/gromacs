#ifndef __SE_FGG_H__
#define __SE_FGG_H__

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "pme-internal.h"
#include "emmintrin.h"
#include "immintrin.h"
#include "x86intrin.h"

// --------------------------------------------------------------------------
static int half(int p)
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
static void
SE_fp_set_zero(real* x, const int N)
{
  memset(x,0.0,N*sizeof(real));
}


// -----------------------------------------------------------------------------
// unpacking params
inline static void
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
inline static void
SE_FGG_FCN_params(SE_FGG_params* params, const SE_opt* opt, int N)
{
  params->N = N;
  params->P = (int) opt->P;
  params->P_half=half( opt->P );
  params->c = opt->c;
  params->d = pow(params->c/PI,1.5);
  params->h = opt->box[0]/opt->M;
  params->a = -FGG_INF;
  params->eta = opt->eta;
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


// ------------------------------------------------------------------------------
void static umr(real *H, int np,char* str)
{
  int d1;real rsum=0;
  for(d1=0;d1<np;d1++)
    rsum += H[np]*H[np];
  printf("%s real %10.6f %g \n",str,rsum,rsum);
}

// ------------------------------------------------------------------------------
void static sumc(t_complex* H, int np,char* str)
{
  int d1;real rsum=0,isum=0;
  for(d1=0;d1<np;d1++)
    if(!isnan(H[np].re) && !isnan(H[np].im))
      {
	rsum += pow(H[np].re,2);
	isum += pow(H[np].im,2);
      }
  printf("%s complex %10.6f %g \n",str,rsum,rsum);
  printf("%s complex %10.6f %g \n",str,isum,isum);
}

// -----------------------------------------------------------------------------
static void
SE_FGG_base_gaussian(real* zs, const SE_FGG_params* params)
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

  //#ifdef _OPENMP
  //#pragma omp for private(i,j,k) schedule(static)// work-share over OpenMP threads here
  //#endif
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
	      zs[idx++] = d*exp(-c*ijkh2);
	    }
	}
    }
}



// -----------------------------------------------------------------------------
static int
fgg_expansion_all(const real x[3], const real q,
                  const SE_FGG_params* params,
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
  const      real h = params->h;
  const      real c = params->c;

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
          idx = (int) round(x[j]/h);
          t0[j] = x[j]-h*idx;
        }
    }
  else
    {
      for(j=0; j<3; j++)
        {
          idx = (int) floor(x[j]/h);
          t0[j] = x[j]-h*idx;
        }
    }

  // compute third factor 
  real z3 = exp(-c*(t0[0]*t0[0] + t0[1]*t0[1] + t0[2]*t0[2]) )*q;

  // compute second factor by induction
  real z_base0 = exp(2*c*h*t0[0]);
  real z_base1 = exp(2*c*h*t0[1]);
  real z_base2 = exp(2*c*h*t0[2]);

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

real static sesum(real *f, int n, int e1, int e2, int dim, char* str)
{
  real s1=0,s2=0;
  int i,j,k;
  if(dim==1)
    for (i=e1;i<e2;i++){
      s1+=f[i];
      s2+=f[i]*f[i];
    }
  else if(dim==3)
    for (i=e1;i<e2;i++)
      for (j=e1;j<e2;j++)
	for (k=e1;k<e2;k++){
	  s1+=f[i*n*n+j*n+k];
	  s2+=f[i*n*n+j*n+k]*f[i*n*n+j*n+k];
	}
  printf("%s %g %s^2 %f\n",str,s1,str,s2);
  return s2;

}


#endif //__SE_FGG_H__
