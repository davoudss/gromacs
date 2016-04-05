#ifndef _SE_GRID_AVX_256_SINGLE_
#define _SE_GRID_AVX_256_SINGLE_

/* SE AVX 256 single gridding */
#ifndef GMX_DOUBLE
#ifdef GMX_X86_AVX_256

#include "se.h"

// -----------------------------------------------------------------------------
void SE_grid_split_AVX(real* grid, real* q,
		       splinedata_t *spline, 
		       const SE_FGG_params* params)
{
  // unpack parameters
  const int     N = params->N;
  float*        H = (float*) grid; // pointer to grid does NOT alias
  const float* zs = (float*) spline->zs;
  const float* zx = (float*) spline->theta[0];
  const float* zy = (float*) spline->theta[1];
  const float* zz = (float*) spline->theta[2];
    
  const int p = params->P;
  float qn;
  int idx0, idx_zs, idx_zz, i, j, k, n, nn;
  const int incrj = params->npdims[2]-p; // middle increment
  const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

  __m256 rH0, rZZ0, rZS0, rC;

  for(n=0; n<N; n++){
    nn = spline->ind[n];
    qn = q[nn];
    idx0 = spline->idx[n];
    idx_zs = 0;

    if(idx0%8 == 0){ // H[idx0] is 32-aligned
      for(i = 0; i<p; i++){
	for(j = 0; j<p; j++){
	  rC = _mm256_set1_ps( qn*zx[p*n+i]*zy[p*n+j] );
	  idx_zz=p*n;
	  for(k = 0; k<p; k+=8){
	    rH0  = _mm256_load_ps( H+idx0     );
	    rZZ0 = _mm256_load_ps( zz + idx_zz     );
	    rZS0 = _mm256_load_ps( zs + idx_zs    );

	    rZZ0 = _mm256_mul_ps(rZZ0,rC);
	    rZZ0 = _mm256_mul_ps(rZZ0,rZS0);
	    rH0  = _mm256_add_ps(rH0,rZZ0);

	    _mm256_store_ps( H+idx0    , rH0 );

	    idx0  +=8;
	    idx_zs+=8; 
	    idx_zz+=8;
	  }
	  idx0 += incrj; 
	}
	idx0 += incri; 
      }
    }
    else{ // H[idx0] is 16-aligned, preventing nice vectorization
      for(i = 0; i<p; i++){
	for(j = 0; j<p; j++){
	  rC = _mm256_set1_ps( qn*zx[p*n+i]*zy[p*n+j] );
	  idx_zz=p*n;
	  for(k = 0; k<p; k+=8){
	    rH0  = _mm256_loadu_ps( H+idx0 );

	    rZZ0 = _mm256_load_ps( zz + idx_zz );
	    rZS0 = _mm256_load_ps( zs + idx_zs );

	    rZZ0 = _mm256_mul_ps(rZZ0,rC);
	    rZZ0 = _mm256_mul_ps(rZZ0,rZS0);

	    rH0  = _mm256_add_ps(rH0,rZZ0);
	    _mm256_storeu_ps( H+idx0, rH0 );

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
void 
SE_grid_split_AVX_dispatch(real* grid, real* q, 
			   splinedata_t *spline,
			   const SE_FGG_params* params)
{
  const int p = params->P;
  const int incrj = params->dims[2]; // middle increment
  const int incri = params->npdims[2]*(params->dims[1]);// outer increment

  // if P, or either increments are not divisible by 4, fall back on vanilla
  if( isnot_div_by_8(p) || isnot_div_by_8(incri) || isnot_div_by_8(incrj) || (p%8)!=0)
    {
      __DISPATCHER_MSG("[FGG GRID AVX] AVX Abort (PARAMS)\n");
      SE_grid_split_SSE_dispatch(grid, q, spline, params);
      return;
    }
    
  // otherwise the preconditions for AVX codes are satisfied. 
  else if(p==8){
    // specific for p=8
    __DISPATCHER_MSG("[FGG GRID AVX] P=8\n");
    SE_grid_split_AVX(grid, q, spline, params); 
  }
}

#endif //GMX_X86_AVX_256
#endif //not GMX_DOUBLE
#endif // _SE_GRID_AVX_256_SINGLE_
