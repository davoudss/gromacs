

/*   THIS CODE IS EDITED BY DAVOUD SAFFAR  */

// =====================================================================
// GROMACS headers
// ======================================================================

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "vec.h"
#include "gmxcomplex.h"
#include "smalloc.h"
#include "futil.h"
#include "gmx_fatal.h"
#include "physics.h"
#include "coulomb.h"
#include <fftw3.h>

#define __FFT
//#define __FFTW


//#ifdef __FFT
#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREAD_MPI
#include "tmpi.h"
#endif

#include "smalloc.h"
#include "gmx_parallel_3dfft.h"
#include "gmx_fft.h"
#include "gmxcomplex.h"
#include "gmx_fatal.h"
#include "fft5d.h"

typedef struct gmx_se {
  MPI_Comm     mpi_comm;
  int          nthread;
  real        *fftgridA; 
  int          fftgrid_nx, fftgrid_ny, fftgrid_nz;
  t_complex   *cfftgridA;
  int          cfftgrid_nx, cfftgrid_ny, cfftgrid_nz;
  int          segrid_nx,segrid_ny,segrid_nz;
  
  gmx_parallel_3dfft_t  pfft_setupA;
} t_gmx_se;

typedef t_gmx_se *gmx_se_t;

//#endif

//#define __POTENTIAL__
#define __FORCE__

#define TOL 2e-4

// =====================================================================
// SE headers
// ======================================================================

#include "SE_fgg.h"

#define PI 3.141592653589793

// Malloc with 16-byte alignment (from intrinsics library)
#define SE_FGG_MALLOC(sz) _mm_malloc((sz),16)
#define SE_FGG_FREE(sz) _mm_free((sz))


// =====================================================================
// SE FGG Utility Routines
// ======================================================================

inline int is_odd(int p)
// test if interger is odd or even
{
  return p&1;
}

// Return half of gaussian width:
//    q = (p-1)/2 if p is odd
//    q =  p/2    if p is even
inline int half(int p)
{
  return (is_odd(p) ? (p-1)/2 : p/2);
}


// -----------------------------------------------------------------------------
// vector index mod, e.g. for real x[6], x[7] --(vmod(7,6))--> x[1] 
static inline int vmod(int i, int N)
{
  if(i>=0)
    return i%N;

  int k = -i/N;
  i = i + (k+1)*N;
  return i % N;
}

// -----------------------------------------------------------------------------
static real randnum(real min, real L)
{
  real q = ( (real) rand() )/RAND_MAX;
  return L*q+min;
}

// =============================================================================
// SE FGG Utility routines =====================================================

// -----------------------------------------------------------------------------
// Set array elements to real-precision zero
static void 
SE_fp_set_zero(real* x, const int N)
{
  int i;
  for(i=0; i<N; i++)
    x[i] = 0.0;
}

// -----------------------------------------------------------------------------
static inline int SE_prod3(const int v[3])
{
  return v[0]*v[1]*v[2];
}

// -----------------------------------------------------------------------------
inline real SE_gettime(void) 
{
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return tv.tv_sec + 1e-6*tv.tv_usec;
}

// -----------------------------------------------------------------------------
static void 
SE_FGG_pack_params(SE_FGG_params* params, int N, int M0, int M1, int M2, 
		   int P, real c, real h)
{
  params->N = N;
  params->P = P;
  params->P_half=half(P);
  params->c = c;
  params->d = pow(c/PI,1.5);
  params->h = h;
  params->a = -1;

  params->dims[0] = M0;
  params->dims[1] = M1;
  params->dims[2] = M2;

  params->npdims[0] = M0+P;
  params->npdims[1] = M1+P;
  params->npdims[2] = M2+P;
}

// -----------------------------------------------------------------------------
static void 
SE_FGG_allocate_workspace(SE_FGG_work* work, const SE_FGG_params* params, 
			  int allocate_zs, int allocate_fgg_expa)
{
  const int P=params->P;
  int numel = SE_prod3(params->npdims);
  work->H = SE_FGG_MALLOC(numel*sizeof(real));
  SE_fp_set_zero(work->H, numel);

  if(allocate_zs)
    work->zs = SE_FGG_MALLOC(P*P*P*sizeof(real));
  else
    work->zs = NULL;

  work->free_zs=allocate_zs;
    
  if(allocate_fgg_expa)
    {
      numel = (params->N)*(params->P);
      work->zx = (real*) SE_FGG_MALLOC(numel*sizeof(real));
      work->zy = (real*) SE_FGG_MALLOC(numel*sizeof(real));;
      work->zz = (real*) SE_FGG_MALLOC(numel*sizeof(real));;
      work->idx= (int*) SE_FGG_MALLOC(params->N*sizeof(int));
    }
  else
    {
      work->zx=NULL;
      work->zy=NULL;
      work->zz=NULL;
      work->idx=NULL;
    }
  work->free_fgg_expa=allocate_fgg_expa;
}

// -----------------------------------------------------------------------------
static real* 
SE_FGG_allocate_grid(const SE_FGG_params* params)
{
  int numel = SE_prod3(params->dims);
  real* H_per = SE_FGG_MALLOC(numel*sizeof(real));
  SE_fp_set_zero(H_per, numel);
  return H_per;
}

// -----------------------------------------------------------------------------
static real* 
SE_FGG_allocate_vec(const int Nm)
{
  real* phi = SE_FGG_MALLOC(Nm*sizeof(real));
  SE_fp_set_zero(phi, Nm);
  return phi;
}


// -----------------------------------------------------------------------------
static void 
SE_FGG_free_workspace(SE_FGG_work* work)
{
  SE_FGG_FREE(work->H);

  if(work->free_zs)
    SE_FGG_FREE(work->zs);

  if(work->free_fgg_expa)
    {
      SE_FGG_FREE(work->zx);
      SE_FGG_FREE(work->zy);
      SE_FGG_FREE(work->zz);
      SE_FGG_FREE(work->idx);
    }
}


// -----------------------------------------------------------------------------
static void 
SE_init_system(SE_state* s, const SE_FGG_params* params)
{
  int i;
  const int N=params->N;

  s->x = SE_FGG_MALLOC(3*N*sizeof(real));
  s->q = SE_FGG_MALLOC(  N*sizeof(real));

}

// -----------------------------------------------------------------------------
static void 
SE_free_system(SE_state* s)
{
  SE_FGG_FREE(s->x);
  SE_FGG_FREE(s->q);
}

// -----------------------------------------------------------------------------
// Wrap H to produce periodicity
static void 
SE_FGG_wrap_fcn(real* H_per, 
		const SE_FGG_work* work, const SE_FGG_params* params)
{
  int idx ,i,j,k;
  int widx[3];
  const int p_half = half(params->P);

  // can not openMP here, race to += on H_per beacuse indices wrap around
  for(i=0; i<params->npdims[0]; i++)
    {
      for(j=0; j<params->npdims[1]; j++)
	{
	  for(k=0; k<params->npdims[2]; k++)
	    {
	      widx[0] = vmod(i-p_half,params->dims[0]);
	      widx[1] = vmod(j-p_half,params->dims[1]);
	      widx[2] = vmod(k-p_half,params->dims[2]);

	      idx = __IDX3_RMAJ(widx[2], widx[0], widx[1], 
				params->dims[1], params->dims[2]);
	      H_per[idx] += work->H[ __IDX3_RMAJ(i,j,k,
						 params->npdims[1],
						 params->npdims[2])];

	    }
	}
    }
}

// -----------------------------------------------------------------------------
// Extend periodic function larger box
// INPUT IN FORTRAN/MATLAB-STYLE COLUMN MAJOR LAYOUT!
static void 
SE_FGG_extend_fcn(SE_FGG_work* work, const real* H_per, 
		  const SE_FGG_params* params)
{
  int idx ,i,j,k;
  int widx[3];
  const int p_half = params->P_half;

#ifdef _OPENMP
#pragma omp for private(i,j,k) schedule(static)// work-share over OpenMP threads here
#endif
  for(i=0; i<params->npdims[0]; i++)
    {
      for(j=0; j<params->npdims[1]; j++)
	{
	  for(k=0; k<params->npdims[2]; k++)
	    {
	      widx[0] = vmod(i-p_half,params->dims[0]);
	      widx[1] = vmod(j-p_half,params->dims[1]);
	      widx[2] = vmod(k-p_half,params->dims[2]);

	      idx = __IDX3_RMAJ(widx[2], widx[0], widx[1], 
				params->dims[1], params->dims[2]);
	      work->H[__IDX3_RMAJ(i,j,k,params->npdims[1],params->npdims[2])]
		= H_per[idx];
	    }
	}
    }
}



// =============================================================================
// Core SE FGG routines ========================================================

// -----------------------------------------------------------------------------
static void
SE_FGG_base_gaussian(SE_FGG_work* work, const SE_FGG_params* params)
{
  int idx ,i,j,k;
    
  // unpack parameters
  const int p=params->P;
  const int p_half = half(p);
  const int p_from = (is_odd(p) ? p_half:p_half-1);
  const real h=params->h;
  const real c=params->c;
  const real d=params->d;

#ifdef _OPENMP
#pragma omp for private(i,j,k) schedule(static)// work-share over OpenMP threads here
#endif
  for(i = -p_from; i<=p_half; i++)
    {
      // hoisting this index calculation (more) breaks omp-parallel code
      idx = __IDX3_RMAJ(i+p_from, 0, 0, p, p);
      for(j = -p_from; j<=p_half; j++)
	for(k = -p_from; k<=p_half; k++)
	  {
	    work->zs[idx++] = d*exp(-c*((i*h)*(i*h) + 
					(j*h)*(j*h) + 
					(k*h)*(k*h)));
	  }
    }
}

// -----------------------------------------------------------------------------
static int 
fgg_expansion_3p(const real x[3], const real q,
		 const SE_FGG_params* params,
		 real z2_0[P_MAX], 
		 real z2_1[P_MAX], 
		 real z2_2[P_MAX])
{
  // unpack params
  const int p = params->P;
  const int p_half = params->P_half;
  const real h = params->h;
  const real c=params->c;
    
  real t0[3];
  int idx ,i,j;
  int idx_from[3];

  // compute index range and centering
  if(is_odd(p))
    {
      for(j=0; j<3; j++)
	{
	  idx = (int) round(x[j]/h);
	  idx_from[j] = idx - p_half;
	  t0[j] = x[j]-h*idx;
	}
    }
  else
    {
      for(j=0; j<3; j++)
	{
	  idx = (int) floor(x[j]/h);
	  idx_from[j] = idx - (p_half-1);
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
    }	
  else
    {
      z0 = pow(z_base0,-p_half+1);
      z1 = pow(z_base1,-p_half+1);
      z2 = pow(z_base2,-p_half+1);
    }

  z2_0[0] = z0;
  z2_1[0] = z1;
  z2_2[0] = z2;
  for(i=1; i<p; i++)
    {
      z0 *=z_base0;
      z1 *=z_base1;
      z2 *=z_base2;

      z2_0[i] = z0;
      z2_1[i] = z1;
      z2_2[i] = z2;
    }

  // save some flops by multiplying one vector with z3 factor
  for(i=0; i<p; i++)
    {
      z2_0[i] *= z3;
    }

  return __IDX3_RMAJ(idx_from[0]+p_half, 
		     idx_from[1]+p_half, 
		     idx_from[2]+p_half, 
		     params->npdims[1], params->npdims[2]);
}


// -----------------------------------------------------------------------------
static int 
fgg_expansion_3p_force(const real x[3], const real q,
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
    const real h      = params->h;
    const real c      = params->c;
    
    real t0[3];
    int i,j,idx;
    int idx_from[3], p_from;

    // compute index range and centering
    if(is_odd(p))
    {
	for(j=0; j<3; j++)
	{
	    idx = (int) round(x[j]/h);
	    idx_from[j] = idx - p_half;
	    t0[j] = x[j]-h*idx;
	}
    }
    else
    {
	for(j=0; j<3; j++)
	{
	    idx = (int) floor(x[j]/h);
	    idx_from[j] = idx - (p_half-1);
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

	zf_0[i] = -c*(t0[0]-(p_from+i)*h);
	zf_1[i] = -c*(t0[1]-(p_from+i)*h);
	zf_2[i] = -c*(t0[2]-(p_from+i)*h);
    }
 
    // save some flops by multiplying one vector with z3 factor
    for(i=0; i<p; i++)
    {
	z2_0[i] *= z3;
    }

    return __IDX3_RMAJ(idx_from[0]+p_half, 
		       idx_from[1]+p_half, 
		       idx_from[2]+p_half, 
		       params->npdims[1], params->npdims[2]);
}



// -----------------------------------------------------------------------------
static void 
SE_FGG_expand_all(SE_FGG_work* work, 
		  const SE_state* st, 
		  const SE_FGG_params* params)
{
  int n;
  real xn[3] __attribute__((aligned(16)));
  const int N = params->N;
  const int P = params->P;

  for(n=0; n<N; n++)
    {
      // compute index and expansion vectors
      xn[0] = st->x[n]; xn[1] = st->x[n+N]; xn[2] = st->x[n+2*N];
	
      *(work->idx+n) = __FGG_EXPA(xn,1,params, 
				  work->zx+n*P, 
				  work->zy+n*P, 
				  work->zz+n*P);
    }
}

// -----------------------------------------------------------------------------
// vanilla grid gather
static void 
SE_FGG_int(real* phi,  
	   const SE_FGG_work* work, 
	   const SE_state* st, 
	   const SE_FGG_params* params)
{
  real z2_0[P_MAX] __attribute__((aligned(16)));
  real z2_1[P_MAX] __attribute__((aligned(16)));
  real z2_2[P_MAX] __attribute__((aligned(16)));

  // unpack params
  const real* H = work->H;
  const real* zs = work->zs;
  const int p = params->P;
  const int N = params->N;
  const real h=params->h;

  real xm[3];
  int i,j,k,m,idx, zidx;
  real phi_m, cij;

  const int incrj = params->npdims[2]-p;
  const int incri = params->npdims[2]*(params->npdims[1]-p);

#ifdef _OPENMP
#pragma omp for private(m) schedule(static)// work-share over OpenMP threads here
#endif
  for(m=0; m<N; m++)
    {
      xm[0] = st->x[m]; xm[1] = st->x[m+N]; xm[2] = st->x[m+2*N];

      idx = __FGG_EXPA(xm, 1, params, z2_0, z2_1, z2_2);
	
      phi_m = 0;
      zidx = 0;

      for(i = 0; i<p; i++)
	{
	  for(j = 0; j<p; j++)
	    {
	      cij = z2_0[i]*z2_1[j];
	      for(k = 0; k<p; k++)
		{
		  phi_m += H[idx]*zs[zidx]*z2_2[k]*cij;
		  idx++; zidx++;
		}
	      idx += incrj;
	    }
	  idx += incri;
	}
      phi[m] = (h*h*h)*phi_m;
    }
}

// -----------------------------------------------------------------------------
// vanilla grid gather to calculate forces
static void 
SE_FGG_int_force(real* force,  
		 const SE_FGG_work* work, 
		 const SE_state* st, 
		 const SE_FGG_params* params)
{
    real z2_0[P_MAX] __attribute__((aligned(16)));
    real z2_1[P_MAX] __attribute__((aligned(16)));
    real z2_2[P_MAX] __attribute__((aligned(16)));

    // to alculate forces
    real zf_0[P_MAX] __attribute__((aligned(16)));
    real zf_1[P_MAX] __attribute__((aligned(16)));
    real zf_2[P_MAX] __attribute__((aligned(16)));

    // unpack params
    const real* H  = work->H;
    const real* zs = work->zs;
    const int   p  = params->P;
    const int   N  = params->N;
    const real  h  = params->h;

    real xm[3], qm;
    int i,j,k,idx, zidx,m;
    real force_m[3], cij;
    
    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

#ifdef _OPENMP
#pragma omp for private(m) schedule(static)// work-share over OpenMP threads here
#endif
    for(m=0; m<N; m++)
    {
      xm[0] = st->x[m]; xm[1] = st->x[m+N]; xm[2] = st->x[m+2*N]; 
      qm = st->q[m];

      idx = __FGG_EXPA_FORCE(xm, qm, params, z2_0, z2_1, z2_2,zf_0,zf_1,zf_2);
      
      force_m[0] = 0; force_m[1] = 0; force_m[2] = 0;
      zidx = 0;
      
      for(i = 0; i<p; i++)
	{
	  for(j = 0; j<p; j++)
	    {
	      cij = z2_0[i]*z2_1[j];
	      for(k = 0; k<p; k++)
		{
		  force_m[0] += H[idx]*zs[zidx]*z2_2[k]*cij*zf_0[i];
		  force_m[1] += H[idx]*zs[zidx]*z2_2[k]*cij*zf_1[j];
		  force_m[2] += H[idx]*zs[zidx]*z2_2[k]*cij*zf_2[k];
		  
		  idx++; zidx++;
		}
	      idx += incrj;
	    }
	  idx += incri;
	}
      force[m    ] = (h*h*h)*force_m[0];
      force[m+  N] = (h*h*h)*force_m[1];
      force[m+2*N] = (h*h*h)*force_m[2];
    }
}


// -----------------------------------------------------------------------------
void SE_FGG_int_split_SSE_dispatch(real*  phi,  
				   const SE_FGG_work* work, 
				   const SE_FGG_params* params)
{
  const int p = params->P;
  const int incrj = params->dims[2]; // middle increment
  const int incri = params->npdims[2]*(params->dims[1]);// outer increment

#if 0
  // THIS BYPASSES THE FAST SSE KERNELS.
  // 
  // THEY ARE PLATFORM DEPENDENT, AND THUS MAY NOT WORK OUT OF THE BOX.
  // REMOVE THIS BLOCK ONLY IF YOU ARE FAMILIAR WITH BASIC DEBUGGING, 
  // THE BASICS OF SSE INTRINSICS, AND ARE WILLING TO UNDERSTAND WHERE
  // THE (DATA ALIGNMENT) PRECONDITIONS OF SSE INSTRUCTIONS MAY BREAK
  // IN THE SSE CODE BELOW.
  __DISPATCHER_MSG("[FGG INT SSE] SSE Disabled\n");
  SE_FGG_int_split(phi, work, params);
  return;
#endif

  // if P is odd, or if either increment is odd, fall back on vanilla
  if( is_odd(p) || is_odd(incri) || is_odd(incrj) )
    {
      __DISPATCHER_MSG("[FGG INT SSE] SSE Abort (PARAMS)\n");
      SE_FGG_int_split(phi, work, params);
      return;
    }

#if 0
  // If the work arrays zs or zX are misaligned, fall back on vanilla.
  // These arrays are dynamically allocated, so getting this alignment
  // is really the compilers job! Once you trust it, remove this 
  // check, because the long integer modulus operation is not fast.
  if( ( (unsigned long) work->zs)%16 != 0 || 
      ( (unsigned long) work->zx)%16 != 0 || 
      ( (unsigned long) work->zy)%16 != 0 ||
      ( (unsigned long) work->zz)%16 != 0 )
    {
      __DISPATCHER_MSG("[FGG INT SSE] SSE Abort (DATA)\n");
      SE_FGG_int_split(phi, work, params);
      return;
    }
#endif

  // otherwise the preconditions for SSE codes are satisfied. 
  if(p==8)
    {
      // specific for p=8
      __DISPATCHER_MSG("[FGG INT SSE] P=8\n");
      SE_FGG_int_split_SSE_P8(phi, work, params);
    } 
  else if(p==16)
    {
      // specific for p=16
      __DISPATCHER_MSG("[FGG INT SSE] P=16\n");
      SE_FGG_int_split_SSE_P16(phi, work, params); 
    }
  else if(p%8==0)
    {
      // for p divisible by 8
      __DISPATCHER_MSG("[FGG INT SSE] P unroll 8\n");
      SE_FGG_int_split_SSE_u8(phi, work, params); 
    }
  else
    {
      // vanilla SSE code (any even p)
      __DISPATCHER_MSG("[FGG INT SSE] Vanilla\n");
      SE_FGG_int_split_SSE(phi, work, params);
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split(real* phi,  
		      const SE_FGG_work* work, 
		      const SE_FGG_params* params)
{
  // unpack params
  const real* H = work->H;
  const real* zs = work->zs;
  const real* zx = work->zx;
  const real* zy = work->zy;
  const real* zz = work->zz;

  const int p = params->P;
  const int N = params->N;
  const real h=params->h;

  int i,j,k,m,idx,idx_zs,idx_zz;
  real phi_m, cij;

  const int incrj = params->npdims[2]-p;
  const int incri = params->npdims[2]*(params->npdims[1]-p);

#ifdef _OPENMP
#pragma omp for private(m)// work-share over OpenMP threads here
#endif
  for(m=0; m<N; m++)
    {
      idx = work->idx[m];	
      phi_m = 0;
      idx_zs = 0;

      for(i = 0; i<p; i++)
	{
	  for(j = 0; j<p; j++)
	    {
	      cij = zx[m*p+i]*zy[m*p+j];
	      idx_zz=m*p;
	      for(k = 0; k<p; k++)
		{
		  phi_m += H[idx]*zs[idx_zs]*zz[idx_zz]*cij;
		  idx++; idx_zs++; idx_zz++;
		}
	      idx += incrj;
	    }
	  idx += incri;
	}
      phi[m] = (h*h*h)*phi_m;
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split_SSE(real* phi,  
			  const SE_FGG_work* work, 
			  const SE_FGG_params* params)
{
  // unpack params
  const real* H = work->H;
  const real* zs = work->zs;
  const real* zx = work->zx;
  const real* zy = work->zy;
  const real* zz = work->zz;

  const int p = params->P;
  const int N = params->N;
  const real h=params->h;

  int i,j,k,m,idx,idx_zs,idx_zz;
  real s[2] __attribute__((aligned(16)));

  __m128d rH0, rZZ0, rZS0, rC, rP;

  const int incrj = params->npdims[2]-p;
  const int incri = params->npdims[2]*(params->npdims[1]-p);

  for(m=0; m<N; m++)
    {
      idx = work->idx[m];	
      idx_zs = 0;
      rP=_mm_setzero_pd();

      if(idx%2==0) // H[idx] is 16-aligned so vectorization simple
	{
	  for(i = 0; i<p; i++)
	    {
	      for(j = 0; j<p; j++)
		{
		  rC = _mm_set1_pd( zx[m*p+i]*zy[m*p+j]);
		  idx_zz=m*p;
		  for(k = 0; k<p; k+=2)
		    {
		      rH0  = _mm_load_pd( H+idx );
		      rZZ0 = _mm_load_pd( zz + idx_zz);
		      rZS0 = _mm_load_pd( zs + idx_zs);
		      rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
			
		      idx+=2; 
		      idx_zs+=2; 
		      idx_zz+=2;
		    }
		  idx += incrj;
		}
	      idx += incri;
	    }
	}
      else // H[idx] not 16-aligned, so use non-aligned loads
	{
	  for(i = 0; i<p; i++)
	    {
	      for(j = 0; j<p; j++)
		{
		  rC = _mm_set1_pd( zx[m*p+i]*zy[m*p+j]);
		  idx_zz=m*p;
		  for(k = 0; k<p; k+=2)
		    {
		      rH0  = _mm_loadu_pd( H+idx );
		      rZZ0 = _mm_load_pd( zz + idx_zz);
		      rZS0 = _mm_load_pd( zs + idx_zs);
		      rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
			
		      idx+=2; 
		      idx_zs+=2; 
		      idx_zz+=2;
		    }
		  idx += incrj;
		}
	      idx += incri;
	    }

	}
      _mm_store_pd(s,rP);
      phi[m] = (h*h*h)*(s[0]+s[1]);
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split_SSE_P8(real*   phi,  
			     const SE_FGG_work* work, 
			     const SE_FGG_params* params)
{
  // unpack params
  const real*   H = work->H;
  const real*   zs = work->zs;
  const real*   zx = work->zx;
  const real*   zy = work->zy;
  const real*   zz = work->zz;
  
  /* ASSUME P=8 const int p = params->P; */
  const int N = params->N;
  const real h=params->h;
  
  int i,j,m,idx,idx_zs;
  real s[2] __attribute__((aligned(16)));
  
  // hold entire zz vector
  __m128d rZZ0, rZZ1, rZZ2, rZZ3; 
  __m128d rC, rP;
  __m128d rH0, rH1, rH2, rH3; 
  __m128d rZS0, rZS1, rZS2, rZS3;
  
  const int incrj = params->npdims[2]-8;
  const int incri = params->npdims[2]*(params->npdims[1]-8);

  for(m=0; m<N; m++)
    {
      idx = work->idx[m];
      idx_zs = 0;
      rP=_mm_setzero_pd();
      
      /* hoist load of ZZ vector */
      rZZ0 = _mm_load_pd(zz + m*8     );
      rZZ1 = _mm_load_pd(zz + m*8 + 2 );
      rZZ2 = _mm_load_pd(zz + m*8 + 4 );
      rZZ3 = _mm_load_pd(zz + m*8 + 6 );
      
      if(idx%2==0) // H[idx] is 16-aligned so vectorization simple
	{
	  for(i = 0; i<8; i++)
	    {
	      for(j = 0; j<8; j++)
		{
		  rC = _mm_set1_pd( zx[m*8+i]*zy[m*8+j]);
		  
		  rH0  = _mm_load_pd( H+idx    );
		  rH1  = _mm_load_pd( H+idx + 2);
		  rH2  = _mm_load_pd( H+idx + 4);
		  rH3  = _mm_load_pd( H+idx + 6);
		  
		  rZS0 = _mm_load_pd( zs + idx_zs    );
		  rZS1 = _mm_load_pd( zs + idx_zs + 2);
		  rZS2 = _mm_load_pd( zs + idx_zs + 4);
		  rZS3 = _mm_load_pd( zs + idx_zs + 6);
		  
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1)));
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2)));
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3)));
		  
		  idx_zs +=8;
		  idx += incrj + 8;
		}
	      idx += incri;
	    }
	}
      else // H[idx] not 16-aligned, so use non-aligned loads
	{
	  for(i = 0; i<8; i++)
	    {
	      for(j = 0; j<8; j++)
		{
		  rC = _mm_set1_pd( zx[m*8+i]*zy[m*8+j]);

		  rH0  = _mm_loadu_pd( H+idx    );
		  rH1  = _mm_loadu_pd( H+idx + 2);
		  rH2  = _mm_loadu_pd( H+idx + 4);
		  rH3  = _mm_loadu_pd( H+idx + 6);

		  rZS0 = _mm_load_pd( zs + idx_zs    );
		  rZS1 = _mm_load_pd( zs + idx_zs + 2);
		  rZS2 = _mm_load_pd( zs + idx_zs + 4);
		  rZS3 = _mm_load_pd( zs + idx_zs + 6);
		    
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1)));
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2)));
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3)));

		  idx_zs +=8;
		  idx += incrj + 8;
		}
	      idx += incri;
	    }
	}
      _mm_store_pd(s,rP);
      phi[m] = (h*h*h)*(s[0]+s[1]);
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split_SSE_P16(real*   phi,  
			      const SE_FGG_work* work, 
			      const SE_FGG_params* params)
{
  // unpack params
  const real*   H = work->H;
  const real*   zs = work->zs;
  const real*   zx = work->zx;
  const real*   zy = work->zy;
  const real*   zz = work->zz;

  /* ASSUME P=16 const int p = params->P; */
  const int N = params->N;
  const real h=params->h;

  int i,j,m,idx,idx_zs;
  real s[2] __attribute__((aligned(16)));

  // hold entire zz vector
  __m128d rZZ0, rZZ1, rZZ2, rZZ3, rZZ4, rZZ5, rZZ6, rZZ7; 
  __m128d rC, rP;
  __m128d rH0, rZS0;

  const int incrj = params->npdims[2]-16;
  const int incri = params->npdims[2]*(params->npdims[1]-16);

  for(m=0; m<N; m++)
    {
      idx = work->idx[m];
      _mm_prefetch( (void*) (H+idx), _MM_HINT_T0);

      idx_zs = 0;
      _mm_prefetch( (void*) zs, _MM_HINT_T0);

      rP=_mm_setzero_pd();

      /* hoist load of ZZ vector */
      rZZ0 = _mm_load_pd(zz + m*16     );
      rZZ1 = _mm_load_pd(zz + m*16 + 2 );
      rZZ2 = _mm_load_pd(zz + m*16 + 4 );
      rZZ3 = _mm_load_pd(zz + m*16 + 6 );
      rZZ4 = _mm_load_pd(zz + m*16 + 8 );
      rZZ5 = _mm_load_pd(zz + m*16 + 10);
      rZZ6 = _mm_load_pd(zz + m*16 + 12);
      rZZ7 = _mm_load_pd(zz + m*16 + 14);
      
      if(idx%2==0) // H[idx] is 16-aligned so vectorization simple
	{
	  for(i = 0; i<16; i++)
	    {
	      for(j = 0; j<16; j++)
		{
		  rC = _mm_set1_pd( zx[m*16+i]*zy[m*16+j]);
		  
		  /* 0 */ 
		  rH0  = _mm_load_pd( H+idx );
		  rZS0 = _mm_load_pd( zs + idx_zs);		    
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));

		  /* 1 */ 
		  rH0  = _mm_load_pd( H+idx + 2);
		  rZS0 = _mm_load_pd( zs + idx_zs + 2);		    
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS0)));

		  /* 2 */ 
		  rH0  = _mm_load_pd( H+idx + 4);
		  rZS0 = _mm_load_pd( zs + idx_zs + 4);		    
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS0)));

		  /* 3 */ 
		  rH0  = _mm_load_pd( H+idx + 6);
		  rZS0 = _mm_load_pd( zs + idx_zs + 6);		    
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS0)));

		  /* 4 */ 
		  rH0  = _mm_load_pd( H+idx + 8);
		  rZS0 = _mm_load_pd( zs + idx_zs + 8);		    
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rC),rZS0)));

		  /* 5 */ 
		  rH0  = _mm_load_pd( H+idx + 10);
		  rZS0 = _mm_load_pd( zs + idx_zs + 10);		    
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rC),rZS0)));

		  /* 6 */ 
		  rH0  = _mm_load_pd( H+idx + 12);
		  rZS0 = _mm_load_pd( zs + idx_zs + 12);		    
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rC),rZS0)));

		  /* 7 */ 
		  rH0  = _mm_load_pd( H+idx + 14);
		  rZS0 = _mm_load_pd( zs + idx_zs + 14);		    
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rC),rZS0)));

		  idx_zs +=16;
		  idx += incrj + 16;
		}
	      idx += incri;
	    }
	}
      else // H[idx] not 16-aligned, so use non-aligned loads
	{
	  for(i = 0; i<16; i++)
	    {
	      for(j = 0; j<16; j++)
		{
		  rC = _mm_set1_pd( zx[m*16+i]*zy[m*16+j]);

		  /* 0 */ 
		  rH0  = _mm_loadu_pd( H+idx );
		  rZS0 = _mm_load_pd( zs + idx_zs);		    
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));

		  /* 1 */ 
		  rH0  = _mm_loadu_pd( H+idx + 2);
		  rZS0 = _mm_load_pd( zs + idx_zs + 2);		    
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS0)));

		  /* 2 */ 
		  rH0  = _mm_loadu_pd( H+idx + 4);
		  rZS0 = _mm_load_pd( zs + idx_zs + 4);		    
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS0)));

		  /* 3 */ 
		  rH0  = _mm_loadu_pd( H+idx + 6);
		  rZS0 = _mm_load_pd( zs + idx_zs + 6);		    
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS0)));

		  /* 4 */ 
		  rH0  = _mm_loadu_pd( H+idx + 8);
		  rZS0 = _mm_load_pd( zs + idx_zs + 8);		    
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rC),rZS0)));

		  /* 5 */ 
		  rH0  = _mm_loadu_pd( H+idx + 10);
		  rZS0 = _mm_load_pd( zs + idx_zs + 10);		    
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rC),rZS0)));

		  /* 6 */ 
		  rH0  = _mm_loadu_pd( H+idx + 12);
		  rZS0 = _mm_load_pd( zs + idx_zs + 12);		    
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rC),rZS0)));

		  /* 7 */ 
		  rH0  = _mm_loadu_pd( H+idx + 14);
		  rZS0 = _mm_load_pd( zs + idx_zs + 14);		    
		  rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rC),rZS0)));

		  idx_zs +=16;
		  idx += incrj + 16;
		}
	      idx += incri;
	    }
	}
      _mm_store_pd(s,rP);
      phi[m] = (h*h*h)*(s[0]+s[1]);
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split_SSE_u8(real*   phi,  
			     const SE_FGG_work* work, 
			     const SE_FGG_params* params)
{
  // unpack params
  const real*   H = work->H;
  const real*   zs = work->zs;
  const real*   zx = work->zx;
  const real*   zy = work->zy;
  const real*   zz = work->zz;

  const int p = params->P;
  const int N = params->N;
  const real h=params->h;

  int i,j,k,m,idx,idx_zs,idx_zz;
  real s[2] __attribute__((aligned(16)));

  __m128d rH0, rZZ0, rZS0, rC, rP;
  __m128d rH1, rZZ1, rZS1;
  __m128d rH2, rZZ2, rZS2;
  __m128d rH3, rZZ3, rZS3;

  const int incrj = params->npdims[2]-p;
  const int incri = params->npdims[2]*(params->npdims[1]-p);

  for(m=0; m<N; m++)
    {
      idx = work->idx[m];
      _mm_prefetch( (void*) (H+idx), _MM_HINT_T0);
	
      idx_zs = 0;
      _mm_prefetch( (void*) zs, _MM_HINT_T0);

      rP=_mm_setzero_pd();

      if(idx%2==0) // H[idx] is 16-aligned so vectorization simple
	{
	  for(i = 0; i<p; i++)
	    {
	      for(j = 0; j<p; j++)
		{
		  rC = _mm_set1_pd( zx[m*p+i]*zy[m*p+j] );
		  idx_zz=m*p;

		  for(k = 0; k<p; k+=8)
		    {
		      rH0  = _mm_load_pd( H+idx    );
		      rH1  = _mm_load_pd( H+idx + 2);
		      rH2  = _mm_load_pd( H+idx + 4);
		      rH3  = _mm_load_pd( H+idx + 6);

		      rZZ0 = _mm_load_pd( zz + idx_zz    );
		      rZZ1 = _mm_load_pd( zz + idx_zz + 2);
		      rZZ2 = _mm_load_pd( zz + idx_zz + 4);
		      rZZ3 = _mm_load_pd( zz + idx_zz + 6);

		      rZS0 = _mm_load_pd( zs + idx_zs    );
		      rZS1 = _mm_load_pd( zs + idx_zs + 2);
		      rZS2 = _mm_load_pd( zs + idx_zs + 4);
		      rZS3 = _mm_load_pd( zs + idx_zs + 6);

		      rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
		      rP = _mm_add_pd(rP,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1)));
		      rP = _mm_add_pd(rP,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2)));
		      rP = _mm_add_pd(rP,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3)));
			
		      idx+=8; 
		      idx_zs+=8; 
		      idx_zz+=8;
		    }
		  idx += incrj;
		}
	      idx += incri;
	    }
	}
      else // H[idx] not 16-aligned, so use non-aligned load from H
	{
	  for(i = 0; i<p; i++)
	    {
	      for(j = 0; j<p; j++)
		{
		  rC = _mm_set1_pd( zx[m*p+i]*zy[m*p+j] );
		  idx_zz=m*p;

		  for(k = 0; k<p; k+=8)
		    {
		      rH0  = _mm_loadu_pd( H+idx    );
		      rH1  = _mm_loadu_pd( H+idx + 2);
		      rH2  = _mm_loadu_pd( H+idx + 4);
		      rH3  = _mm_loadu_pd( H+idx + 6);

		      rZZ0 = _mm_load_pd( zz + idx_zz    );
		      rZZ1 = _mm_load_pd( zz + idx_zz + 2);
		      rZZ2 = _mm_load_pd( zz + idx_zz + 4);
		      rZZ3 = _mm_load_pd( zz + idx_zz + 6);

		      rZS0 = _mm_load_pd( zs + idx_zs    );
		      rZS1 = _mm_load_pd( zs + idx_zs + 2);
		      rZS2 = _mm_load_pd( zs + idx_zs + 4);
		      rZS3 = _mm_load_pd( zs + idx_zs + 6);

		      rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
		      rP = _mm_add_pd(rP,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1)));
		      rP = _mm_add_pd(rP,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2)));
		      rP = _mm_add_pd(rP,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3)));
			
		      idx+=8; 
		      idx_zs+=8; 
		      idx_zz+=8;
		    }
		  idx += incrj;
		}
	      idx += incri;
	    }
	}

      // done accumulating
      _mm_store_pd(s,rP);
      phi[m] = (h*h*h)*(s[0]+s[1]);
    }
}

// -----------------------------------------------------------------------------
static void 
SE_FGG_grid(SE_FGG_work* work, const SE_state* st, 
	    const SE_FGG_params* params)
{
  // vectors for FGG expansions
  real zx0[P_MAX] __attribute__((aligned(16)));
  real zy0[P_MAX] __attribute__((aligned(16)));
  real zz0[P_MAX] __attribute__((aligned(16)));

  // unpack parameters
  const int N=params->N;
  real*   H = work->H; // pointer to grid does NOT alias
  const real*   zs = work->zs;
  const int p = params->P;
    
  real cij0;
  real xn[3];
  real qn;
  int idx0, zidx, i,j,k,n;
  const int incrj = params->npdims[2]-p; // middle increment
  const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

#ifdef _OPENMP
#pragma omp for private(n)// work-share over OpenMP threads here
#endif
  for(n=0; n<N; n++)
    {
      // compute index and expansion vectors
      xn[0] = st->x[n]; xn[1] = st->x[n+N]; xn[2] = st->x[n+2*N];
      qn = st->q[n];
      idx0 = __FGG_EXPA(xn, qn, params, zx0, zy0, zz0);
      zidx = 0;
      for(i = 0; i<p; i++)
	{
	  for(j = 0; j<p; j++)
	    {
	      cij0 = zx0[i]*zy0[j];
	      for(k = 0; k<p; k++)
		{
		  H[idx0] += zs[zidx]*zz0[k]*cij0;
		  idx0++; 
		  zidx++;
		}
	      idx0 += incrj; 
	    }
	  idx0 += incri; 
	}
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split_SSE_dispatch(SE_FGG_work* work, const SE_state* st, 
				    const SE_FGG_params* params)
{
  const int p = params->P;
  const int incrj = params->dims[2]; // middle increment
  const int incri = params->npdims[2]*(params->dims[1]);// outer increment

#if 0
  // THIS BYPASSES THE FAST SSE KERNELS.
  // 
  // THEY ARE PLATFORM DEPENDENT, AND THUS MAY NOT WORK OUT OF THE BOX.
  // REMOVE THIS BLOCK ONLY IF YOU ARE FAMILIAR WITH BASIC DEBUGGING, 
  // THE BASICS OF SSE INTRINSICS, AND ARE WILLING TO UNDERSTAND WHERE
  // THE (DATA ALIGNMENT) PRECONDITIONS OF SSE INSTRUCTIONS MAY BREAK
  // IN THE SSE CODE BELOW.
  __DISPATCHER_MSG("[FGG GRID SSE] SSE Disabled\n");
  SE_FGG_grid_split(work, st, params);
  return;
#endif

  // if P is odd, or if either increment is odd, fall back on vanilla
  if( is_odd(p) || is_odd(incri) || is_odd(incrj) )
    {
      __DISPATCHER_MSG("[FGG GRID SSE] SSE Abort (PARAMS)\n");
      SE_FGG_grid_split(work, st, params);
      return;
    }

#if 0
  // If the work arrays zs or zX are misaligned, fall back on vanilla.
  // These arrays are dynamically allocated, so getting this alignment
  // is really the compilers job! Once you trust it, remove this 
  // check, because the long integer modulus operation is not fast.
  if( ( (unsigned long) work->zs)%16 != 0 || 
      ( (unsigned long) work->zx)%16 != 0 || 
      ( (unsigned long) work->zy)%16 != 0 ||
      ( (unsigned long) work->zz)%16 != 0 )
    {
      __DISPATCHER_MSG("[FGG GRID SSE] SSE Abort (DATA)\n");
      SE_FGG_grid_split(work, st, params);
      return;
    }
#endif

  // otherwise the preconditions for SSE codes are satisfied. 
  if(p==16)
    {
      // specific for p=16
      __DISPATCHER_MSG("[FGG GRID SSE] P=16\n");
      SE_FGG_grid_split_SSE_P16(work, st, params); 
    }
  else if(p%8==0)
    {
      // specific for p divisible by 8
      __DISPATCHER_MSG("[FGG GRID SSE] P unroll 8\n");
      SE_FGG_grid_split_SSE_u8(work, st, params); 
    }
  else
    {
      // vanilla SSE code (any even p)
      __DISPATCHER_MSG("[FGG GRID SSE] Vanilla\n");
      SE_FGG_grid_split_SSE(work, st, params);
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split(SE_FGG_work* work, const SE_state* st, 
		       const SE_FGG_params* params)
{
  // unpack parameters
  const int N=params->N;
  real*   H = work->H; // pointer to grid does NOT alias
  const real*   zs = work->zs;
  const real*   zx = work->zx;
  const real*   zy = work->zy;
  const real*   zz = work->zz;
    
  const int p = params->P;

  real cij0;
  real qn;
  int idx0, zidx, idxzz, i, j, k,n;
  const int incrj = params->npdims[2]-p; // middle increment
  const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

#ifdef _OPENMP
#pragma omp for private(n)// work-share over OpenMP threads here
#endif
  for(n=0; n<N; n++)
    {
      qn = st->q[n];
      idx0 = work->idx[n];

      // inline vanilla loop
      zidx = 0;
      for(i = 0; i<p; i++)
	{
	  for(j = 0; j<p; j++)
	    {
	      cij0 = qn*zx[p*n+i]*zy[p*n+j];
	      idxzz=p*n;
	      for(k = 0; k<p; k++)
		{
		  H[idx0] += zs[zidx]*zz[idxzz]*cij0;
		  idx0++; zidx++; idxzz++;
		}
	      idx0 += incrj; 
	    }
	  idx0 += incri; 
	}
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split_SSE_P16(SE_FGG_work* work, const SE_state* st, 
			       const SE_FGG_params* params)
{
  // unpack parameters
  const int N=params->N;
  real*   H = work->H; // pointer to grid does NOT alias
  const real*   zs = work->zs;
  const real*   zx = work->zx;
  const real*   zy = work->zy;
  const real*   zz = work->zz;

  real qn;
  int idx, idx_zs, i,j,n;
  const int incrj = params->npdims[2]-16; // middle increment
  const int incri = params->npdims[2]*(params->npdims[1]-16);// outer increment

  __m128d rZZ0, rZZ1, rZZ2, rZZ3, rZZ4, rZZ5, rZZ6, rZZ7; 
  __m128d rH0, rH1, rH2, rH3;
  __m128d rC, rZS0;

  for(n=0; n<N; n++)
    {
      qn = st->q[n];

      idx = work->idx[n];
      _mm_prefetch( (void*) (H+idx), _MM_HINT_T0);

      idx_zs = 0;
      _mm_prefetch( (void*) zs, _MM_HINT_T0);

      rZZ0 = _mm_load_pd(zz + n*16     );
      rZZ1 = _mm_load_pd(zz + n*16 + 2 );
      rZZ2 = _mm_load_pd(zz + n*16 + 4 );
      rZZ3 = _mm_load_pd(zz + n*16 + 6 );
      rZZ4 = _mm_load_pd(zz + n*16 + 8 );
      rZZ5 = _mm_load_pd(zz + n*16 + 10);
      rZZ6 = _mm_load_pd(zz + n*16 + 12);
      rZZ7 = _mm_load_pd(zz + n*16 + 14);

      if(idx%2 == 0) // H[idx0] is 16-aligned
	{
	  for(i = 0; i<16; i++)
	    {
	      for(j = 0; j<16; j++)
		{
		  rC = _mm_set1_pd( qn*zx[16*n+i]*zy[16*n+j] );

		  /* 0 - 3 */ 
		  rH0  = _mm_load_pd( H+idx    );
		  rH1  = _mm_load_pd( H+idx + 2);
		  rH2  = _mm_load_pd( H+idx + 4);
		  rH3  = _mm_load_pd( H+idx + 6);

		  rZS0 = _mm_load_pd( zs + idx_zs);
		  rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0));

		  rZS0 = _mm_load_pd( zs + idx_zs + 2);                   
		  rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS0));

		  rZS0 = _mm_load_pd( zs + idx_zs + 4);                   
		  rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS0));

		  rZS0 = _mm_load_pd( zs + idx_zs + 6);                   
		  rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS0));

		  _mm_store_pd(H + idx, rH0);
		  _mm_store_pd(H + idx + 2, rH1);
		  _mm_store_pd(H + idx + 4, rH2);
		  _mm_store_pd(H + idx + 6, rH3);

		  /* 4 - 7*/ 
		  rH0  = _mm_load_pd( H+idx + 8 );
		  rH1  = _mm_load_pd( H+idx + 10);
		  rH2  = _mm_load_pd( H+idx + 12);
		  rH3  = _mm_load_pd( H+idx + 14);

		  rZS0 = _mm_load_pd( zs + idx_zs + 8);
		  rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rC),rZS0));

		  rZS0 = _mm_load_pd( zs + idx_zs + 10);                   
		  rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ5,rC),rZS0));

		  rZS0 = _mm_load_pd( zs + idx_zs + 12);                   
		  rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ6,rC),rZS0));

		  rZS0 = _mm_load_pd( zs + idx_zs + 14);                   
		  rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ7,rC),rZS0));

		  _mm_store_pd(H + idx + 8 , rH0);
		  _mm_store_pd(H + idx + 10, rH1);
		  _mm_store_pd(H + idx + 12, rH2);
		  _mm_store_pd(H + idx + 14, rH3);

		  idx += incrj + 16;
		  idx_zs += 16;
		}
	      idx += incri;
	    }
	}
      else // H[idx0] is 8-aligned, preventing nice vectorization
	{
	  for(i = 0; i<16; i++)
	    {
	      for(j = 0; j<16; j++)
		{
		  rC = _mm_set1_pd( qn*zx[16*n+i]*zy[16*n+j] );

		  /* 0 - 3 */ 
		  rH0  = _mm_loadu_pd( H+idx    );
		  rH1  = _mm_loadu_pd( H+idx + 2);
		  rH2  = _mm_loadu_pd( H+idx + 4);
		  rH3  = _mm_loadu_pd( H+idx + 6);

		  // if zs does not have 16-byte alignment, this will core.
		  // PLATFORM AND COMPILER DEPENDENT (FIXME)
		  rZS0 = _mm_load_pd( zs + idx_zs);
		  rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0));

		  rZS0 = _mm_load_pd( zs + idx_zs + 2);                   
		  rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS0));

		  rZS0 = _mm_load_pd( zs + idx_zs + 4);                   
		  rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS0));

		  rZS0 = _mm_load_pd( zs + idx_zs + 6);                   
		  rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS0));

		  _mm_storeu_pd(H + idx, rH0);
		  _mm_storeu_pd(H + idx + 2, rH1);
		  _mm_storeu_pd(H + idx + 4, rH2);
		  _mm_storeu_pd(H + idx + 6, rH3);

		  /* 4 - 7*/ 
		  rH0  = _mm_loadu_pd( H+idx + 8 );
		  rH1  = _mm_loadu_pd( H+idx + 10);
		  rH2  = _mm_loadu_pd( H+idx + 12);
		  rH3  = _mm_loadu_pd( H+idx + 14);

		  rZS0 = _mm_load_pd( zs + idx_zs + 8);
		  rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rC),rZS0));

		  rZS0 = _mm_load_pd( zs + idx_zs + 10);                   
		  rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ5,rC),rZS0));

		  rZS0 = _mm_load_pd( zs + idx_zs + 12);                   
		  rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ6,rC),rZS0));

		  rZS0 = _mm_load_pd( zs + idx_zs + 14);                   
		  rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ7,rC),rZS0));

		  _mm_storeu_pd(H + idx + 8 , rH0);
		  _mm_storeu_pd(H + idx + 10, rH1);
		  _mm_storeu_pd(H + idx + 12, rH2);
		  _mm_storeu_pd(H + idx + 14, rH3);

		  idx += incrj + 16;
		  idx_zs += 16;
		}
	      idx += incri;
	    }
	}
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split_SSE_u8(SE_FGG_work* work, const SE_state* st, 
			      const SE_FGG_params* params)
{
  // unpack parameters
  const int N=params->N;
  real*   H = work->H; // pointer to grid does NOT alias
  const real*   zs = work->zs;
  const real*   zx = work->zx;
  const real*   zy = work->zy;
  const real*   zz = work->zz;
    
  const int p = params->P;

  real qn;
  int idx0, idx_zs, idx_zz, i,j,k,n;
  const int incrj = params->npdims[2]-p; // middle increment
  const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

  __m128d rH0, rZZ0, rZS0, rC;
  __m128d rH1, rZZ1, rZS1;
  __m128d rH2, rZZ2, rZS2;
  __m128d rH3, rZZ3, rZS3;

  for(n=0; n<N; n++)
    {
      qn = st->q[n];

      idx0 = work->idx[n];
      _mm_prefetch( (void*) (H+idx0), _MM_HINT_T0);

      idx_zs = 0;
      _mm_prefetch( (void*) zs, _MM_HINT_T0);

      if(idx0%2 == 0) // H[idx0] is 16-aligned
	{
	  for(i = 0; i<p; i++)
	    {
	      for(j = 0; j<p; j++)
		{
		  rC = _mm_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
		  idx_zz=p*n;

		  for(k = 0; k<p; k+=8)
		    {
		      rH0  = _mm_load_pd( H+idx0     );
		      rH1  = _mm_load_pd( H+idx0 + 2 );
		      rH2  = _mm_load_pd( H+idx0 + 4 );
		      rH3  = _mm_load_pd( H+idx0 + 6 );

		      rZZ0 = _mm_load_pd( zz + idx_zz     );
		      rZZ1 = _mm_load_pd( zz + idx_zz + 2 );
		      rZZ2 = _mm_load_pd( zz + idx_zz + 4 );
		      rZZ3 = _mm_load_pd( zz + idx_zz + 6 );

		      rZS0 = _mm_load_pd( zs + idx_zs    );
		      rZS1 = _mm_load_pd( zs + idx_zs + 2);
		      rZS2 = _mm_load_pd( zs + idx_zs + 4);
		      rZS3 = _mm_load_pd( zs + idx_zs + 6);

		      rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0));
		      rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1));
		      rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2));
		      rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3));
			
		      _mm_store_pd( H+idx0    , rH0 );
		      _mm_store_pd( H+idx0 + 2, rH1 );
		      _mm_store_pd( H+idx0 + 4, rH2 );
		      _mm_store_pd( H+idx0 + 6, rH3 );

		      idx0  +=8;
		      idx_zs+=8; 
		      idx_zz+=8;
		    }
		  idx0 += incrj;
		}
	      idx0 += incri;
	    }
	}
      else // H[idx0] is 8-aligned, preventing nice vectorization
	{
	  for(i = 0; i<p; i++)
	    {
	      for(j = 0; j<p; j++)
	    	{
		  rC = _mm_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
		  idx_zz=p*n;
		    
		  for(k = 0; k<p; k+=8)
	    	    {
		      rH0  = _mm_loadu_pd( H+idx0     );
		      rH1  = _mm_loadu_pd( H+idx0 + 2 );
		      rH2  = _mm_loadu_pd( H+idx0 + 4 );
		      rH3  = _mm_loadu_pd( H+idx0 + 6 );

		      rZZ0 = _mm_load_pd( zz + idx_zz     );
		      rZZ1 = _mm_load_pd( zz + idx_zz + 2 );
		      rZZ2 = _mm_load_pd( zz + idx_zz + 4 );
		      rZZ3 = _mm_load_pd( zz + idx_zz + 6 );

		      rZS0 = _mm_load_pd( zs + idx_zs    );
		      rZS1 = _mm_load_pd( zs + idx_zs + 2);
		      rZS2 = _mm_load_pd( zs + idx_zs + 4);
		      rZS3 = _mm_load_pd( zs + idx_zs + 6);

		      rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0));
		      rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1));
		      rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2));
		      rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3));

		      _mm_storeu_pd( H+idx0    , rH0 );
		      _mm_storeu_pd( H+idx0 + 2, rH1 );
		      _mm_storeu_pd( H+idx0 + 4, rH2 );
		      _mm_storeu_pd( H+idx0 + 6, rH3 );

		      idx0  +=8;
		      idx_zs+=8;
		      idx_zz+=8;
	    	    }
		  idx0 += incrj;
	    	}
	      idx0 += incri;
	    }
	}
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split_SSE(SE_FGG_work* work, const SE_state* st, 
			   const SE_FGG_params* params)
{
  // unpack parameters
  const int N=params->N;
  real*   H = work->H; // pointer to grid does NOT alias
  const real*   zs = work->zs;
  const real*   zx = work->zx;
  const real*   zy = work->zy;
  const real*   zz = work->zz;
    
  const int p = params->P;

  real qn;
  int idx0, idx_zs, idx_zz, i,j,k,n;
  const int incrj = params->npdims[2]-p; // middle increment
  const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

  __m128d rH0, rZZ0, rZS0, rC;

  for(n=0; n<N; n++)
    {
      qn = st->q[n];
      idx0 = work->idx[n];
      idx_zs = 0;

      if(idx0%2 == 0) // H[idx0] is 16-aligned
	{
	  for(i = 0; i<p; i++)
	    {
	      for(j = 0; j<p; j++)
		{
		  rC = _mm_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
		  idx_zz=p*n;
		  for(k = 0; k<p; k+=2)
		    {
		      rH0  = _mm_load_pd( H+idx0     );
		      rZZ0 = _mm_load_pd( zz + idx_zz     );
		      rZS0 = _mm_load_pd( zs + idx_zs    );

		      rZZ0 = _mm_mul_pd(rZZ0,rC);
		      rZZ0 = _mm_mul_pd(rZZ0,rZS0);
		      rH0  = _mm_add_pd(rH0,rZZ0);

		      _mm_store_pd( H+idx0    , rH0 );

		      idx0  +=2;
		      idx_zs+=2; 
		      idx_zz+=2;
		    }
		  idx0 += incrj; 
		}
	      idx0 += incri; 
	    }
	}
      else // H[idx0] is 8-aligned, preventing nice vectorization
	{
	  for(i = 0; i<p; i++)
	    {
	      for(j = 0; j<p; j++)
	    	{
		  rC = _mm_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
		  idx_zz=p*n;
		  for(k = 0; k<p; k+=2)
	    	    {
		      rH0  = _mm_loadu_pd( H+idx0 );
		      rZZ0 = _mm_load_pd( zz + idx_zz );
		      rZS0 = _mm_load_pd( zs + idx_zs );
		      rZZ0 = _mm_mul_pd(rZZ0,rC);
		      rZZ0 = _mm_mul_pd(rZZ0,rZS0);
		      rH0  = _mm_add_pd(rH0,rZZ0);
		      _mm_storeu_pd( H+idx0, rH0 );

		      idx0  +=2;
		      idx_zs+=2;
		      idx_zz+=2;
	    	    }
		  idx0 += incrj;
	    	}
	      idx0 += incri;
	    }
	}
    }
}

// =============================================================================
// Reordering ==================================================================

// -----------------------------------------------------------------------------
static
int compare_idx_pair(const void* a, const void* b)
{
  int temp = ( * ((idx_reorder_t*) a) ).idx_on_grid - 
    ( * ((idx_reorder_t*) b) ).idx_on_grid;
  if (temp > 0)
    return 1;
  else if (temp < 0)
    return -1;
  else
    return 0;
}

// -----------------------------------------------------------------------------
static
SE_state* SE_clone_state(const SE_state* s, const SE_FGG_params* params)
{
  const int N = params->N;

  SE_state* t = (SE_state*) SE_FGG_MALLOC( sizeof(SE_state) );
  t->x = (real*) SE_FGG_MALLOC(3*N*sizeof(real));
  t->q = (real*) SE_FGG_MALLOC(  N*sizeof(real));
    
  memcpy(t->x, s->x, 3*N*sizeof(real));
  memcpy(t->q, s->q,   N*sizeof(real));

  return t;
}

// -----------------------------------------------------------------------------
static void 
SE_FGG_reorder_system(SE_state* s, 
		      const SE_FGG_work* work, 
		      const SE_FGG_params* params)
{
  int i;
  const int N = params->N;

  /* pairing of grid-index and index in array */
  idx_reorder_t* order = (idx_reorder_t*) SE_FGG_MALLOC(N*sizeof(idx_reorder_t));
    
  /* fill index pairing */
  for(i=0;i<N;i++)
    {
      order[i].idx_on_grid=work->idx[i];
      order[i].idx_in_array=i;
    }

  /* sort the temporary */
  qsort(order, N, sizeof(idx_reorder_t), compare_idx_pair);

  /* copy system */
  SE_state* s_copy = SE_clone_state(s, params);

  /* permute system */ 
  int j;
  for(i=0; i<N; i++)
    {
      j = order[i].idx_in_array;

      s->x[i]     = s_copy->x[j  ];
      s->x[i+N]   = s_copy->x[j+N];
      s->x[i+2*N] = s_copy->x[j+2*N];
      s->q[i]     = s_copy->q[j];
    }

  /* deallocations */
  SE_FGG_FREE(order);
  SE_free_system(s_copy); // free contents of s_copy
  SE_FGG_FREE(s_copy);           // s_copy itself is malloc'd
}




// =====================================================================
// SE General Utility Routines
// ======================================================================

// -----------------------------------------------------------------------------
// unpacking params
inline static void 
parse_params(SE_opt* opt, real xi)
{
  //  int      k_inf = 28; // It follows a formula that can be used
  //  int         M0 = 2*k_inf+1;
  //  real      m0 = 0.9*sqrt(PI*(real)(opt->P));  // try also 1.71
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
// create k_space evctors
inline static void
k_vec(int M, real* box,real *k1,real *k2, real *k3)
{
  int i,MM;
  real factor1 = 2.0*PI/box[0];
  real factor2 = 2.0*PI/box[1];
  real factor3 = 2.0*PI/box[2];
  real iter = 0.;
  

  if((M%2)==0){
    MM = M/2;
    for (i=0;i<=MM-1;i++)
      {
	k1[i] = factor1*iter;
	k2[i] = factor2*iter;
	k3[i] = factor3*iter;
	iter += 1.;
      }
		
    iter = (real) 0.-MM;
    for (i=0;i<=MM-1;i++)
      {
	k1[MM+i] = factor1*iter;
	k2[MM+i] = factor2*iter;
	k3[MM+i] = factor3*iter;
	iter += 1.;
      }
  }
  else {
    MM = (M-1)/2;
    iter = 0.;
    for (i=0;i<=MM;i++)
      {
	k1[i] = factor1*iter;
	k2[i] = factor2*iter;
	k3[i] = factor3*iter;
	iter += 1.;
      }
		
    iter = -MM;
    for (i=0;i<=MM-1;i++)
      {
	k1[MM+1+i] = factor1*iter;
	k2[MM+1+i] = factor2*iter;
	k3[MM+1+i] = factor3*iter;
	iter += 1.;
      }
  }

}



// -----------------------------------------------------------------------------
// creating k square
inline static void
k_square(real *k1, real* k2, real* k3, // input vectors
	 real* K2, 		       // output vector
	 int n1, int n2, int n3)       // dimension
{
  int i,j,k;
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif	
  for(i=0;i<n1;i++)
    for(j=0;j<n2;j++)
      for(k=0;k<n3;k++)
	K2[i*n3*n2+j*n2+k] = k1[i]*k1[i]+k2[j]*k2[j]+k3[k]*k3[k];

}

// -----------------------------------------------------------------------------
// packing SE parameters
static void scaling(real *K2, real scalar, real *Z,
	     int n1, int n2, int n3)
{
  int i,j,k;
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k) schedule(static)
#endif	
  for(i=0;i<n1;i++)
    for(j=0;j<n2;j++)
      for(k=0;k<n3;k++)
	Z[i*n3*n2+j*n2+k] = exp(scalar*K2[i*n3*n2+j*n2+k])/K2[i*n3*n2+j*n2+k];

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

  params->dims[0] = opt->M;
  params->dims[1] = opt->M;
  params->dims[2] = opt->M;

  params->npdims[0] = params->dims[0]+params->P;
  params->npdims[1] = params->dims[1]+params->P;
  params->npdims[2] = params->dims[2]+params->P;

}

// -----------------------------------------------------------------------------
// calling gridding
static void
SE_fg_grid(real* x, real* q, int N, SE_opt opt, real* H_out)
{
  // pack parameters
  SE_FGG_params params;
  SE_FGG_FCN_params(&params, &opt, N);

  // scratch arrays
  SE_FGG_work work;
  SE_FGG_allocate_workspace(&work, &params,true,false);

  // initialize the output with zeros
  SE_fp_set_zero(H_out,SE_prod3(params.dims));

  // coordinates and charges
  const SE_state st = {.x = x, .q = q};
 
  if(VERBOSE)
    printf("[SE FG(G)] N=%d, P=%d\n",N,params.P);

  // now do the work
  SE_FGG_base_gaussian(&work, &params);
  SE_FGG_grid(&work, &st, &params);

  SE_FGG_wrap_fcn(H_out, &work, &params);

  // done  
  SE_FGG_free_workspace(&work);
}


// -----------------------------------------------------------------------------
// interopolation and integration
static void
SE_fgg_int(real* x, real* q, real* H, int N, SE_opt opt, real* phi_force_out)
{
  // pack parameters
  SE_FGG_params params;
  SE_FGG_FCN_params(&params, &opt, N);

  // scratch arrays
  SE_FGG_work work;
  SE_FGG_allocate_workspace(&work, &params,true,false);

  // coordinates and charges
#ifdef __POTENTIAL__
  const SE_state st = {.x = x, .q = NULL};
#endif

#ifdef __FORCE__
  const SE_state st = {.x = x, .q = q};
#endif
 
  if(VERBOSE)
    printf("[SE FG(G)] N=%d, P=%d\n",N,params.P);

  // now do the work
  SE_FGG_base_gaussian(&work, &params);
  SE_FGG_extend_fcn(&work, H, &params);

  // Integration
#ifdef __POTENTIAL__
  SE_FGG_int(phi_force_out, &work, &st, &params);
#endif

#ifdef __FORCE__
  SE_FGG_int_force(phi_force_out, &work, &st, &params);
#endif

  // done  
  SE_FGG_free_workspace(&work);
}


// -----------------------------------------------------------------------------
// 3Dfft using fftw3 complex to complex
static void 
do_fft_c2c_forward_3d(fft_complex* in, fft_complex* out ,
		      int n1, int n2, int n3)
{
  if(sizeof(fft_complex)==2*sizeof(double)){
    fftw_plan p;
    typedef double comp[2];
    p = fftw_plan_dft_3d(n1, n2, n3, (comp*) in,(comp*) out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute_dft(p,(comp*) in,(comp*) out);
    fftw_destroy_plan(p);
  }
  else{
    fftwf_plan p;
    typedef float comp[2];
    p = fftwf_plan_dft_3d(n1, n2, n3, (comp*) in, (comp*) out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute_dft(p,(comp*) in,(comp*) out);
    fftwf_destroy_plan(p);
  }
}


// -----------------------------------------------------------------------------
// 3Dfft using fftw3 complex to complex
static void 
do_fft_c2c_backward_3d(fft_complex* in, fft_complex* out ,
		       int n1, int n2, int n3)
{
  if(sizeof(fft_complex)==2*sizeof(double)){
    fftw_plan p;
    typedef double comp[2];
    p = fftw_plan_dft_3d(n1, n2, n3, (comp*) in, (comp*) out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute_dft(p,(comp*) in,(comp*) out);
    fftw_destroy_plan(p);
  }
  else{
    fftwf_plan p;
    typedef float comp[2];
    p = fftwf_plan_dft_3d(n1, n2, n3,(comp*) in, (comp*) out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute_dft(p,(comp*) in,(comp*) out);
    fftwf_destroy_plan(p);
  }
}


// -----------------------------------------------------------------------------
// products sr(scalar to real). flag = 1 gives a.*b and flag = -1 gives a./b
void
product_sr(real* a, real scalar, real *c, int n1, int n2, int n3,int flag)
{
  int i;

  switch (flag)
    {
    case 1:

#ifdef _OPENMP
#pragma omp parallel for private(i) shared(c)
#endif
      for(i=0;i<n1*n2*n3;i++){
	c[i] = a[i]*scalar;
      }
      break;
    case -1:
#ifdef _OPENMP
#pragma omp parallel for private(i) shared(c)
#endif
      for(i=0;i<n1*n2*n3;i++){
	c[i] = a[i]/scalar;
      }
    }

}


// -----------------------------------------------------------------------------
// products rc (real to complex) equivalent to .* in MATLAB. 
// flag = 1 gives a.*b and flag = -1 gives a./b
void
product_rc(fft_complex* a, real* b, fft_complex *c, int n1,int n2, int n3,int flag)
{
  int i;

  switch (flag)
    {
    case 1: 
#ifdef _OPENMP
#pragma omp parallel for private(i) shared(c)
#endif
      for(i=0;i<n1*n2*n3;i++){
	c[i][0] = a[i][0]*b[i];
	c[i][1] = a[i][1]*b[i];
      }
      break;
    case -1:
#ifdef _OPENMP
#pragma omp parallel for private(i) shared(c)
#endif
      for(i=0;i<n1*n2*n3;i++){
	c[i][0] = a[i][0]/b[i];
	c[i][1] = a[i][1]/b[i];
      }
    }

}

// -----------------------------------------------------------------------------
// products rc (real to complex) equivalent to .* in MATLAB. 
// flag = 1 gives a.*b and flag = -1 gives a./b
void
se_product_rc(t_complex* a, real* b, t_complex *c, int n1,int n2, int n3,int flag)
{
  int i;

  switch (flag)
    {
    case 1: 
#ifdef _OPENMP
#pragma omp parallel for private(i) shared(c)
#endif
      for(i=0;i<n1*n2*n3;i++){
	c[i].re = a[i].re*b[i];
	c[i].im = a[i].im*b[i];
      }
      break;
    case -1:
#ifdef _OPENMP
#pragma omp parallel for private(i) shared(c)
#endif
      for(i=0;i<n1*n2*n3;i++){
	c[i].re = a[i].re/b[i];
	c[i].im = a[i].im/b[i];
      }
    }

}

// -----------------------------------------------------------------------------
// printing results GROMAX rvec 3d array while multiplying by c
void print_rvec(char* str, rvec a[], int n,real c,int verb)
{
  if(verb==0)
    return;
  printf("%s:\n",str);
  int i;
  for (i=0;i<n;i++)
    printf("(%.4f,%.4f,%.4f)\n",a[i][0]*c,a[i][1]*c,a[i][2]*c);
  printf("\n");

}


// -----------------------------------------------------------------------------
// printing results GROMAX real 1d array while multiplying by c
void print_real(char* str, real a[], int n,real c,int verb)
{
  if(verb==0)
    return;
  printf("%s:\n",str);
  int i;
  for (i=0;i<n;i++)
    printf("%.4f\n",a[i]*c);
  printf("\n");

}


// -----------------------------------------------------------------------------
// printing results real 1d array while multiplying by c
void print_r1d(char* str, real *a, int n, real c,int verb)
{
  if(verb == 0)
    return;
  printf("%s:\n",str);
  int i;
  for (i=0;i<n;i++)
    printf("%.4f\n",a[i]*c);
  printf("\n");
}

// -----------------------------------------------------------------------------
// printing results complex 1d array while multiplying by c
void print_c1d(char* str, fft_complex *a, int n, real c,int verb)
{
  if (verb == 0)
    return;
  printf("%s:\n",str);
  int i;
  for (i=0;i<n;i++)
    printf("(%.4f,%.4f)\t",a[i][0]*c,a[i][1]*c);
  printf("\n");
}

// -----------------------------------------------------------------------------
// printing results real 3d array while multiplying by c
void print_r3d(char* str, real *a, int n, real c,int verb)
{
  if(verb == 0)
    return;
  printf("%s:\n",str);
  int i;
  for (i=0;i<n;i++)
    printf("(%3.4f,%3.4f,%3.4f)\n",a[i]*c,a[i+n]*c,a[i+2*n]*c);
  printf("\n");
}

// -----------------------------------------------------------------------------
// printing results complex 3d array while multiplying by c
void print_c3d(char* str, fft_complex *a, int n1, int n2, int n3, real c,int verb)
{
  if (verb == 0)
    return;
  printf("%s:\n",str);
  int i,j,k;
  for (i=0;i<n1;i++){
    for(j=0;j<n2;j++){
      for(k=0;k<n3;k++)
	printf("(%3.4f,%3.4f)\t",a[i*n3*n2+j*n2+k][0]*c,a[i*n3*n2+j*n2+k][1]*c);
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");
}

// -----------------------------------------------------------------------------
// copy real vec into complex vec

void
copy_r2c(real* in, fft_complex* out,int n)
{
  int i;
  for (i=0;i<n;i++){
    out[i][0] = in[i];
    out[i][1] = 0.;
  }
}


// -----------------------------------------------------------------------------
// copy complex vec into real vec

void
copy_c2r(fft_complex* in, real* out,int n)
{
  int i;
  for (i=0;i<n;i++)
    out[i] = in[i][0];

}

// -----------------------------------------------------------------------------
// max3(a,b,c)

inline int
max3(int a, int b, int c )
{
  int t;
  t = (a>b) ? a : b;
  t = (t>c) ? t : c;
  return t;
}


// ----------------------------------------------------------------------------
// lambertw function similar to MATLAB
/* written K M Briggs Keith dot Briggs at bt dot com 97 May 21.  
   Revised KMB 97 Nov 20; 98 Feb 11, Nov 24, Dec 28; 99 Jan 13; 00 Feb 23; 01 Apr 09

   Computes Lambert W function, principal branch.
   See LambertW1.c for -1 branch.

   Returned value W(z) satisfies W(z)*exp(W(z))=z
   test data...
   W(1)= 0.5671432904097838730
   W(2)= 0.8526055020137254914
   W(20)=2.2050032780240599705
   To solve (a+b*R)*exp(-c*R)-d=0 for R, use
   R=-(b*W(-exp(-a*c/b)/b*d*c)+a*c)/b/c
*/

inline static real
LambertW(const double z) {
  int i; 
  const double eps=4.0e-16, em1=0.3678794411714423215955237701614608; 
  double p,e,t,w;
  if (z<-em1 || isinf(z) || isnan(z)) { 
    fprintf(stderr,"LambertW: bad argument %g, exiting.\n",z); exit(1); 
  }
  if (0.0==z) return 0.0;
  if (z<-em1+1e-4) { // series near -em1 in sqrt(q)
    double q=z+em1,r=sqrt(q),q2=q*q,q3=q2*q;
    return 
      -1.0
      +2.331643981597124203363536062168*r
      -1.812187885639363490240191647568*q
      +1.936631114492359755363277457668*r*q
      -2.353551201881614516821543561516*q2
      +3.066858901050631912893148922704*r*q2
      -4.175335600258177138854984177460*q3
      +5.858023729874774148815053846119*r*q3
      -8.401032217523977370984161688514*q3*q;  // error approx 1e-16
  }
  /* initial approx for iteration... */
  if (z<1.0) { /* series near 0 */
    p=sqrt(2.0*(2.7182818284590452353602874713526625*z+1.0));
    w=-1.0+p*(1.0+p*(-0.333333333333333333333+p*0.152777777777777777777777)); 
  } else 
    w=log(z); /* asymptotic */
  if (z>3.0) w-=log(w); /* useful? */
  for (i=0; i<10; i++) { /* Halley iteration */
    e=exp(w); 
    t=w*e-z;
    p=w+1.0;
    t/=e*p-0.5*(p+1.0)*t/p; 
    w-=t;
    if (fabs(t)<eps*(1.0+fabs(w))) return w; /* rel-abs error */
  }
  /* should never get here */
  fprintf(stderr,"LambertW: No convergence at z=%g, exiting.\n",z); 
  exit(1);
}


// -----------------------------------------------------------------------------
// set M, P, m based on the inputs eps, L, N, xi for starting SE

static real
SE_init_params(int* M,int* P,real* m, real eps, real L,
               real xi, real* q, int N)
{
  int i,iP,M_l,P_l=3;
  double eps_half = eps/2.;
  double m_l, App,c=.95;
  real Q=0.;

  // based on Lindbo & Tornberg error estimate
  for (iP=3;iP<P_MAX;iP++){
    m_l   = c*sqrt(PI*iP);
    App = exp(-iP*iP*PI*PI/(2.*m_l*m_l))+4.*erfc(m_l/sqrt(2.));
    if(App<eps_half){
      P_l = iP;
      break;
    }
  }

  for (i=0;i<N;i++)
    Q += q[i]*q[i];

  // Based on Kolafa & Perram error estimate
  double w = 4./((double) 3.*L*L)*((double) Q/pow(PI*xi*eps_half*eps_half,2./3.));
  double lw = LambertW(w);
  double k_inf = sqrt(3.)*xi*L/(2*PI)*sqrt(lw);

  M_l = max3(ceil(2.*k_inf+1),P_l,ceil(sqrt((double) P_l/PI)*xi*L/c));

  *m = m_l;
  *M = M_l;
  *P = P_l;
  return Q;
}



// copu se grids to fftgrids to perform fft
static void
copy_segrid_to_fftgrid(gmx_parallel_3dfft_t pfft_setup, 
		       real* segrid, real* fftgrid,
		       int nx, int ny, int nz)
{
  ivec local_fft_ndata,local_fft_offset,local_fft_size;
  ivec local_se_size;
  int  seidx,ix,iy,iz,fftidx;
  
  gmx_parallel_3dfft_real_limits(pfft_setup,
				 local_fft_ndata,
				 local_fft_offset,
				 local_fft_size);
  local_se_size[0] = nx;
  local_se_size[1] = ny;
  local_se_size[2] = nz;
  
  for (ix = 0; ix < local_fft_ndata[XX]; ix++){
    for (iy = 0; iy < local_fft_ndata[YY]; iy++){
      for (iz = 0; iz < local_fft_ndata[ZZ]; iz++){
	seidx  = ix*(local_se_size[YY]*local_se_size[ZZ])+iy*(local_se_size[ZZ])+iz;
	fftidx = ix*(local_fft_size[YY]*local_fft_size[ZZ])+iy*(local_fft_size[ZZ])+iz;
	fftgrid[fftidx] = segrid[seidx];
      }
    }
  }
}

static void
copy_fftgrid_to_segrid(gmx_parallel_3dfft_t pfft_setup, 
		       const real *fftgrid, real *segrid,
		       int nthread, int thread,
		       int nx, int ny, int nz)
{
    ivec          local_fft_ndata, local_fft_offset, local_fft_size;
    ivec          local_se_size;
    int           ixy0, ixy1, ixy, ix, iy, iz;
    int           seidx, fftidx;

    /* Dimensions should be identical for A/B grid, so we just use A here */
    gmx_parallel_3dfft_real_limits(pfft_setup,
                                   local_fft_ndata,
                                   local_fft_offset,
                                   local_fft_size);

    local_se_size[0] = nx;
    local_se_size[1] = ny;
    local_se_size[2] = nz;
    
    ixy0 = ((thread  )*local_fft_ndata[XX]*local_fft_ndata[YY])/nthread;
    ixy1 = ((thread+1)*local_fft_ndata[XX]*local_fft_ndata[YY])/nthread;

    for (ixy = ixy0; ixy < ixy1; ixy++)
    {
        ix = ixy/local_fft_ndata[YY];
        iy = ixy - ix*local_fft_ndata[YY];

        seidx = (ix*local_se_size[YY] + iy)*local_se_size[ZZ];
        fftidx = (ix*local_fft_size[YY] + iy)*local_fft_size[ZZ];
        for (iz = 0; iz < local_fft_ndata[ZZ]; iz++)
        {
            segrid[seidx+iz] = fftgrid[fftidx+iz];
        }
    }
    
}


// =====================================================================
// SE main call routine
// ======================================================================
#define Z(i,j,k) Z[M*(i*M+j)+k]

int
spectral_ewald(real *x, real *q, SE_opt opt,
	       real xi, real* phi_force)
{

  // parameters and constants
  int M,N,P;
  real m,w,eta,c;

#ifdef __FFT
  // variables to run fft
  gmx_parallel_3dfft_t pfft_setup;
  real                *fftgrid   = NULL;
  t_complex           *cfftgrid  = NULL;
  int                  bReproducible = 1;
  int                  nthreads = 1;
  MPI_Comm             comm[2];
  int                 *minor = NULL, *major = NULL;
  ivec                 ndata;
  gmx_wallcycle_t      wcycle = NULL;
#endif

  // creating and set new parameters
  parse_params(&opt,xi);

  // setting parameters for use
  m = opt.m; w = opt.w; eta = opt.eta; c = opt.c; 
  M = opt.M; N = opt.N; P = opt.P;

  real *H_in = malloc(sizeof(real)*M*M*M);
#ifdef __FFTW
  fft_complex *H_in_comp = malloc(sizeof(fft_complex)*M*M*M); // to copy real to comp for fft
#endif

  // TO GRID FUNCTION
  SE_fg_grid(x,q,N,opt,H_in);

  // TRANSFORM AND SHIFT
#ifdef __FFTW
  copy_r2c(H_in,H_in_comp,M*M*M);
  fft_complex *H_out = malloc(sizeof(fft_complex)*M*M*M);
  do_fft_c2c_forward_3d(H_in_comp,H_out,M,M,M); 
#endif

#ifdef __FFT
  ndata[0] = M;
  ndata[1] = M;
  ndata[2] = M;
  gmx_parallel_3dfft_init(&pfft_setup,ndata,&fftgrid,&cfftgrid,
			  comm,major,minor,bReproducible,nthreads);

  copy_segrid_to_fftgrid(pfft_setup,H_in,fftgrid,M,M,M);
 
  gmx_parallel_3dfft_execute(pfft_setup,GMX_FFT_REAL_TO_COMPLEX,
			     fftgrid,cfftgrid,0,wcycle);
#endif

  //  printf("********** line %d in %s *********** \n",__LINE__,__FILE__);


  // SCALING
  // k-vectors
  real *k1 = malloc(sizeof(real)*M);
  real *k2 = malloc(sizeof(real)*M);
  real *k3 = malloc(sizeof(real)*M);
  k_vec(M,opt.box,k1,k2,k3); 
  // scale
  real *K2 = malloc(sizeof(real)*M*M*M);
  k_square(k1,k2,k3,K2,M,M,M);
  real scalar = -(1.-eta)/(4.*xi*xi);
  real *Z = malloc(sizeof(real)*M*M*M);
  scaling(K2,scalar,Z,M,M,M);
  Z(0,0,0) = 0.;
  int  flag = 1;
#ifdef __FFTW
  product_rc(H_out,Z,H_out,M,M,M,flag);
#endif

#ifdef __FFT
  se_product_rc(cfftgrid,Z,cfftgrid,M,M,M,flag);
#endif

  // INVERSE SHIFT AND INVERSE TRANSFORM
#ifdef __FFTW
  do_fft_c2c_backward_3d(H_out,H_in_comp,M,M,M);
  copy_c2r(H_in_comp,H_in,M*M*M);
#endif

#ifdef __FFT
  gmx_parallel_3dfft_execute(pfft_setup,GMX_FFT_COMPLEX_TO_REAL,
			     cfftgrid,fftgrid,0,wcycle);

  copy_fftgrid_to_segrid(pfft_setup,fftgrid,H_in,1,0,M,M,M);
#endif

  // SPREAD AND INTEGRATE
  SE_fgg_int(x,q,H_in,N,opt,phi_force);

#ifdef __POTENTIAL__ 
  scalar = 4.*PI/(real) (M*M*M);
  product_sr(phi_force,scalar,phi_force,N,1,1,flag);
#endif

#ifdef __FORCE__
  scalar = 4.*PI/(real) (M*M*M);
  product_sr(phi_force,scalar,phi_force,3*N,1,1,flag);
#endif

  __FREE(H_in);
  __FREE(k1); __FREE(k2); __FREE(k3); __FREE(K2); 
  __FREE(Z);

#ifdef __FFTW
  __FREE(H_in_comp);  
  __FREE(H_out);
#endif

  //printf("********** line %d in %s *********** \n",__LINE__,__FILE__);
#ifdef __FFT
  gmx_parallel_3dfft_destroy(pfft_setup); 
  free(fftgrid); free(cfftgrid);
#endif

  return 0;

}

// =====================================================================
// GROMACS routines
// ======================================================================

struct ewald_tab
{
  int        nx, ny, nz, kmax;
  cvec     **eir;
  fft_complex *tab_xy, *tab_qxyz;
};


void init_ewald_tab(ewald_tab_t *et, const t_commrec *cr, const t_inputrec *ir, FILE *fp)
{
  int n;

  snew(*et, 1);
  if (fp)
    {
      fprintf(fp, "Will do ordinary reciprocal space Ewald sum.\n");
    }

  (*et)->nx       = ir->nkx+1;
  (*et)->ny       = ir->nky+1;
  (*et)->nz       = ir->nkz+1;
  (*et)->kmax     = max((*et)->nx, max((*et)->ny, (*et)->nz));
  (*et)->eir      = NULL;
  (*et)->tab_xy   = NULL;
  (*et)->tab_qxyz = NULL;
}



real do_ewald(FILE *log,       gmx_bool bVerbose,
              t_inputrec *ir,
              rvec x[],        rvec f[],
              real chargeA[],  real chargeB[],
              rvec box,
              t_commrec *cr,   int natoms,
              matrix lrvir,    real ewaldcoeff,
              real lambda,     real *dvdlambda,
              ewald_tab_t et)
{
  
  const int N = natoms;
  const real eps = TOL;
 
  int i,M,P;
  real m,L=box[0], xi=ewaldcoeff;
 
  SE_FGG_params params = {.N = N};
  SE_state st;

  // Initialize the system
  SE_init_system(&st, &params);
  
  st.q = chargeA;
   
  for (i=0;i<N;i++){
    st.x[i    ] = x[i][0];
    st.x[i+  N] = x[i][1];
    st.x[i+2*N] = x[i][2];
  }
 
  // Initialize parameters for SE
  real Q = SE_init_params(&M,&P,&m,eps,L,xi,st.q,N);

  SE_opt opt = {.m = m, .box = {L, L, L}, .P = P, .M = M, .N=N};


#ifdef __POTENTIAL__

  // allocate space for potential (output)
  real *phi = malloc(sizeof(real)*N);
  // set it to zero
  SE_fp_set_zero(phi,N);

  spectral_ewald(st.x,st.q,opt,xi,phi);
  
  print_r1d("phi",phi,N,1.,0);
  
  real energy=0.;
  for (i=0;i<N;i++)
    energy += chargeA[i]*phi[i];
  energy = energy/2.;

  SE_free_system(&st);

  if(VERBOSE){
    if(sizeof(double)==sizeof(real))
      printf("\nDouble Precision new 1\n");
    else
      printf("\nSingle Precision new 1\n");
  }
  
  __FREE(phi);

  return energy;

#endif

#ifdef __FORCE__
  // Allocate space for the force
  real *force = malloc(sizeof(real)*3*N);
  
  // Initialize the vector
  SE_fp_set_zero(force,3*N);

  spectral_ewald(st.x,st.q,opt,xi,force);
  
  if(VERBOSE)
    for(i=0;i<N;i++)
      printf("(%f,%f,%f)\n",force[i],force[i+N],force[i+2*N]);
  
    for (i=0;i<N;i++){
      f[i][XX] = force[i];
      f[i][YY] = force[i+N];
      f[i][ZZ] = force[i+2*N];
    }
    
    SE_free_system(&st);
    
    if(VERBOSE){
      if(sizeof(double)==sizeof(real))
	printf("\nDouble Precision new force\n");
      else
	printf("\nSingle Precision new force\n");
    }
    
    __FREE(force);

   
    return 0;
#endif

}
