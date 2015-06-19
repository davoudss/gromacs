/* SE core is written by Dag Lindbo, dag@kth.se
 * General functions are links between SE and GROMACS.
 * Computation of forces is added later to use in GROMACS
 * Davoud Saffar Shamshirgar davoudss@kth.se
 */

#define __SE_FGG_H

#define VERBOSE 1
#ifndef VERBOSE
#define VERBOSE 0
#endif

// System includes
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include "math.h"
#include "math.h"
#include <sys/time.h>
#include "emmintrin.h"
#include <assert.h>
#include <stdbool.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif


// Constants and indexing 
#define PI 3.141592653589793
#define FGG_INF 1.79769e+308

#define __IDX3_RMAJ(II,IJ,IK,N2,N3) ( (II)*(N2)*(N3)+(IJ)*(N3)+(IK) )

#define __FGG_EXPA fgg_expansion
#define __FGG_EXPA_ALL fgg_expansion_all

// Maximal amount of Gaussian support (defined to help the compiler)
#define P_MAX 32

// Specific includes and defines for MEX-file compilation
#define __MALLOC malloc
#define __PRINTF printf
#define __FREE free
#define __IDX __IDX3_RMAJ


// display debug messages in the SSE dispatcher
#if defined VERBOSE
#define __DISPATCHER_MSG(s) __PRINTF(s)
//#define __DISPATCHER_MSG(s) {}
#else
#define __DISPATCHER_MSG(s) {}
#endif


// =====================================================================
// SE typedefs
// ======================================================================

// Temporary, or auxillary arrays
typedef real fft_complex[2];

typedef struct 
{
  real m,c,box[3],h,xi,w,eta;
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
    real box[3];

} SE_FGG_params;

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

enum {
SE_SPREAD, SE_FFT, SE_IFFT, SE_REDIST, SE_SOLVE, SE_GATHER
};

/*
void SE_log(FILE* fp, double t, int logMode)
{
	SE_tot += t;
	switch (logMode)
	{
		case 0:
			fprintf(fp,"Spreading:\t\t %8.4f\n",t);
			return;
		case 1:
			fprintf(fp,"FFT:      \t\t %8.4f\n",t);
			return;
		case 2:
			fprintf(fp,"IFFT:     \t\t %8.4f\n",t);
			return;
		case 3:
			fprintf(fp,"Precomp:  \t\t %8.4f\n",t);
			return;
		case 4:
			fprintf(fp,"SE solve: \t\t %8.4f\n",t);
			return;
		case 5:
			fprintf(fp,"Gathering:\t\t %8.4f\n",t);
			fprintf(fp,"----------------------\n");
			fprintf(fp,"Total:    \t\t %8.4f\n",SE_tot);
			return;
	}
	fclose(fp);		
}
*/

// =============================================================
// SE GENERAL ROUTINES
// =============================================================

// Unpacking params
inline static void parse_params(SE_opt*, real);

// create k_space evctors
inline static void k_vec(int, real*, real*, real*, real*);

// do the scaling
static void scaling(real *, real *, real *, real , real *, int , int, int);


// products sr(scalar to real) rr (real to real) rc (real to complex)
// equivalent to .* in MATLAB. flag = 1 gives a.*b and flag = -1 gives a./b
void product_rc(fft_complex*, real*, fft_complex*, int, int, int, int);
void se_product_rc(t_complex*, real*, t_complex*, int, int, int, int);

// Packing SE parameters
inline static void
SE_FGG_FCN_params(SE_FGG_params*, const SE_opt*, int);

// calling gridding
static void SE_fg_grid(real*, real*, int, SE_opt, real*);

// integration and interpolation
static void SE_fgg_int(real*, real*, real*, int, SE_opt, real*);
// integration and interpolation and calculate forces
static void SE_fgg_int_force(real*, real *, real*, int, SE_opt, real*);

// 3d forward fft using fftw3 complex to complex
static void do_fft_c2c_forward_3d(fft_complex*, fft_complex*, int, int, int);
// 3d backward fft using fftw3 complex to complex
static void do_fft_c2c_backward_3d(fft_complex*, fft_complex*, int, int, int);


void copy_segrid_to_fftwgrid(real*, fft_complex*, int);
void copy_fftwgrid_to_segrid(fft_complex*, real*, int);

static real
SE_init_params(int*, int*, real*, double, real,
               real, real*, int);

inline int max3(int,int,int);
// lambertW function similar to MATLAB
static double lambertW(const double);


// =============================================================
// SE INTERNAL ROUTINES
// =============================================================

// Allocate workspace (malloc)
static void SE_FGG_allocate_workspace(SE_FGG_work*, const SE_FGG_params*, int, int);
static real* SE_FGG_allocate_grid(const SE_FGG_params*);
static real* SE_FGG_allocate_vec(int);

// Free workspace (free)
static void SE_FGG_free_workspace(SE_FGG_work*);

// Particles to grid SINGLE precision
static void SE_FGG_grid(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE_dispatch(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE_P8(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_AVX_dispatch(SE_FGG_work*,const SE_state*,const SE_FGG_params*);
void SE_FGG_grid_split_AVX(SE_FGG_work*, const SE_state*, const SE_FGG_params*);

// Particles to grid DOUBLE precision
void SE_FGG_grid_split_SSE_dispatch_d(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_d(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE_d(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE_u8_d(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE_P16_d(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_AVX_dispatch_d(SE_FGG_work*,const SE_state*,const SE_FGG_params*);
void SE_FGG_grid_split_AVX_P16_d(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_AVX_P8_d(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_AVX_u8_d(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_AVX_d(SE_FGG_work*, const SE_state*, const SE_FGG_params*);


// Compute all FGG expansion vectors
static void SE_FGG_expand_all(SE_FGG_work*, const SE_state*, const SE_FGG_params*);

// Grid to particles SINGLE precision
static void 
SE_FGG_int(rvec*, const SE_FGG_work*, SE_state*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_dispatch(rvec*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split(rvec*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE(rvec*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_P8(rvec*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_AVX_dispatch(rvec*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_AVX(rvec*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);

// Grid to particles DOUBLE precision
void SE_FGG_int_split_SSE_dispatch_d(rvec*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_d(rvec*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_d(rvec*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_u8_d(rvec*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_P8_d(rvec*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_P16_d(rvec*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_AVX_dispatch_d(rvec*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_AVX_P16_d(rvec*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_AVX_P8_d(rvec*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_AVX_u8_d(rvec*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_AVX_d(rvec*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);


// Static Gaussian on P^3-grid
static void SE_FGG_base_gaussian(SE_FGG_work*, const SE_FGG_params*);

// Wrap function to produce periodicity
static void SE_FGG_wrap_fcn(real*, const SE_FGG_work*, const SE_FGG_params*);

// Extend periodic function
static void SE_FGG_extend_fcn(SE_FGG_work*, const real*, const SE_FGG_params*);

// Return product of elements in integer triplet
static int SE_prod3(const int[3]);

// Set first N elements of array to floating-point zero
static void SE_fp_set_zero(real*, int);

// Set first N elements of rvec to floating-point zero
static void SE_fp_set_rvec_zero(rvec*, int);


static void 
PME_SE_FGG_base_gaussian(real*, const SE_FGG_params*);
real sesum(real *f, int n, int e1, int e2, int dim, char*);
