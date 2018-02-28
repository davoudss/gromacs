#ifndef SE_H
#define SE_H

#include <stdbool.h>

/* Datatypes, constants and primary functions for 
 * the Spectral Ewald method */
// Constants and indexing 
#define PI 3.141592653589793
#define FGG_INF 1.79769e+308

#define __IDX3_RMAJ(II,IJ,IK,N2,N3) ( (II)*(N2)*(N3)+(IJ)*(N3)+(IK) )

// Maximal amount of Gaussian support (defined to help the compiler)
#define P_MAX 32

// Specific includes and defines for MEX-file compilation
#define __MALLOC malloc
#define __PRINTF printf
#define __FREE free
#define __IDX __IDX3_RMAJ


// display debug messages in the SSE dispatcher
//#if defined VERBOSE
#define __DISPATCHER_MSG(s) __PRINTF(s)
//#else
//#define __DISPATCHER_MSG(s) {}
//#endif

#if __AVX__
#define MEM_ALIGNED __attribute__((aligned(32)))
#define SE_FGG_MALLOC(sz) _mm_malloc((sz),32)
#elif __SSE4_1__
#define MEM_ALIGNED __attribute__((aligned(16)))
#define SE_FGG_MALLOC(sz) _mm_malloc((sz),16)
#endif
#define SE_FGG_FREE(sz) _mm_free((sz))

#if __SSE4_1__ // any intrinsic instruction should enable precomputation
#define PRECOMP_FGG_EXPA 1
#else
#define PRECOMP_FGG_EXPA 0
#endif

/*
#ifdef _OPENMP
#include <omp.h>
#endif
*/

// =====================================================================
// SE typedefs
// ======================================================================

// Temporary, or auxillary arrays
typedef real fft_complex[2];

typedef struct 
{
  real m,c,box[3],h,xi,w,eta,beta;
  int P,M,N;
} SE_opt;

typedef struct
{
  real* H;
  real* zs;
  
  real* zx;
  real* zy;
  real* zz;
  int* idx;
  
  // extra space to compute the force
  real* zfx;
  real* zfy;
  real* zfz;
  
  int free_zs;
  int free_fgg_expa;
  
} SE_FGG_work;

// FGG parameters
typedef struct
{
    int N;
    int P;
    int P_half;
    int dims[3];
    int npdims[3];
    real c;
    real d;
    real h;
    real a;
    real eta;
    real beta;
    real box[3];

} SE_params;

typedef struct
{
    int idx_on_grid, idx_in_array; 
} idx_reorder_t;

// Particle positions and charges
typedef struct
{
    rvec* x;
    real* q;
    real* phi;

} SE_state;


// test if interger is odd or even
static int is_odd(int p)
{
  return p&1;
}

// test if integer is not divisible by 4 for SINGLE SSE, DOUBLE AVX
inline int isnot_div_by_4(int p)
{
  int divisor = 4, quotient, remainder;

  __asm__ ( "movl   %2, %%edx;"
            "sarl  $31, %%edx;"
            "movl   %2, %%eax;"
            "movl   %3, %%ebx;"
            "idivl      %%ebx;"
	    : "=a" (quotient), "=d" (remainder)
	    : "g"  (p), "g"  (divisor)
	    : "ebx" );
  if(remainder==0)
    return 0;
  else
    return 1;

}
// test if integer is not divisible by 8 for SINGLE AVX
inline int isnot_div_by_8(int p)
{
  int divisor = 8, quotient, remainder;

  __asm__ ( "movl   %2, %%edx;"
            "sarl  $31, %%edx;"
            "movl   %2, %%eax;"
            "movl   %3, %%ebx;"
            "idivl      %%ebx;"
	    : "=a" (quotient), "=d" (remainder)
	    : "g"  (p), "g"  (divisor)
	    : "ebx" );
  if(remainder==0)
    return 0;
  else
    return 1;
}

#endif //SE_H
