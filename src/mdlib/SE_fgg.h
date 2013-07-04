/* SE core is written by Dag Lindbo, dag@kth.se
 * General functions are links between SE and GROMACS.
 * Computation of forces is added later to use in GROMACS
 * Davoud Saffar Shamshirgar davoudss@kth.se
 */

#define __SE_FGG_H

#ifndef VERBOSE
#define VERBOSE 0
#endif

// System includes
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include "math.h"
#include <sys/time.h>
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

//#define __IDX3_CMAJ(II,IJ,IK,N1,N2) ( (II)+(IJ)*(N1)+(IK)*(N1)*(N2) )
#define __IDX3_RMAJ(II,IJ,IK,N2,N3) ( (II)*(N2)*(N3)+(IJ)*(N3)+(IK) )

// Select periodicty: must give -D<...> to compiler
#define __FGG_EXPA fgg_expansion_3p
#define __FGG_EXPA_FORCE fgg_expansion_3p_force

// Maximal amount of Gaussian support (defined to help the compiler)
#define P_MAX 32

// Specific includes and defines for MEX-file compilation

#define __MALLOC malloc
#define __PRINTF printf
#define __FREE free
#define __IDX __IDX3_RMAJ


// display debug messages in the SSE dispatcher
#ifdef VERBOSE
#define __DISPATCHER_MSG(s) __PRINTF(s)
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

} SE_FGG_params;

typedef struct
{
    int idx_on_grid, idx_in_array; 
} idx_reorder_t;

// Particle positions and charges
typedef struct
{
    real* x;
    real* q;

} SE_state;



// =============================================================
// SE GENERAL ROUTINES
// =============================================================

// Unpacking params
inline static void parse_params(SE_opt*, real);

// create k_space evctors
inline static void k_vec(int, real*, real*, real*, real*);

// creating k square
inline static void k_square(real*, real*, real*,real*, int, int, int);

// do the scaling
static void scaling(real *, real , real *, int , int, int);


// products sr(scalar to real) rr (real to real) rc (real to complex)
// equivalent to .* in MATLAB. flag = 1 gives a.*b and flag = -1 gives a./b
void product_sr(real*, real, real*, int, int, int, int);
void product_rr(real*, real*, real*, int, int, int, int);
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
// 3d fft using fftw3 real to complex
//void do_fft_r2c_3d(real*, fft_complex*, int, int, int);
// 3d fft using fftw3 complex to real
//void do_fft_c2r_3d(fft_complex*, real*, int, int, int);
// 3d forward fft using fftw3 complex to complex
static void do_fft_c2c_forward_3d(fft_complex*, fft_complex*, int, int, int);
// 3d backward fft using fftw3 complex to complex
static void do_fft_c2c_backward_3d(fft_complex*, fft_complex*, int, int, int);



// printing results
void print_r1d(char*, real*, int, real,int);
void print_rvec(char*, rvec a[], int, real,int);
void print_real(char*, real a[], int, real,int);
void print_c1d(char*, fft_complex*, int, real,int);
void print_r3d(char*, real*, int, real,int);
void print_c3d(char*, fft_complex*, int, int, int, real,int);

void copy_r2c(real*, fft_complex*, int);
void copy_c2r(fft_complex*, real*, int);

static real
SE_init_params(int*, int*, real*, real, real,
               real, real*, int);

inline int max3(int,int,int);
// lambertW function similar to MATLAB
inline static real lambertW(const real);




// =============================================================
// SE INTERNAL ROUTINES
// =============================================================

// Fill parameter struct
static void SE_FGG_pack_params(SE_FGG_params*, int, int, int, int, int, 
			real, real);

// Allocate workspace (malloc)
static void SE_FGG_allocate_workspace(SE_FGG_work*, const SE_FGG_params*, int, int);
static real* SE_FGG_allocate_grid(const SE_FGG_params*);
static real* SE_FGG_allocate_vec(int);

// Free workspace (free)
static void SE_FGG_free_workspace(SE_FGG_work*);

// Particles to grid
static void SE_FGG_grid(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE_dispatch(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE_u8(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_grid_split_SSE_P16(SE_FGG_work*, const SE_state*, const SE_FGG_params*);

// Compute all FGG expansion vectors
static void SE_FGG_expand_all(SE_FGG_work*, const SE_state*, const SE_FGG_params*);

// Grid to particles
static void 
SE_FGG_int(real*, const SE_FGG_work*, const SE_state*, const SE_FGG_params*);
static void 
SE_FGG_int_force(real*,const SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_dispatch(real*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split(real*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE(real*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_u8(real*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_P8(real*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FGG_int_split_SSE_P16(real*, const SE_FGG_work*, const SE_FGG_params*);

// Static Gaussian on P^3-grid
static void SE_FGG_base_gaussian(SE_FGG_work*, const SE_FGG_params*);

// Wrap function to produce periodicity
static void SE_FGG_wrap_fcn(real*, const SE_FGG_work*, const SE_FGG_params*);

// Extend periodic function
static void SE_FGG_extend_fcn(SE_FGG_work*, const real*, const SE_FGG_params*);

// Randomize positions and charges (malloc)
static void SE_init_system(SE_state*, const SE_FGG_params*);

// Free particles and charges (free)
static void SE_free_system(SE_state*);

// Retrun time in seconds
real SE_gettime(void);

// Return product of elements in integer triplet
static int SE_prod3(const int[3]);

// Set first N elements of array to floating-point zero
static void SE_fp_set_zero(real*, int);

// Reorder particles according to closest grid point
static void SE_FGG_reorder_system(SE_state*, const SE_FGG_work*, const SE_FGG_params*);

