#ifndef _SE_GRID_AVX_256_SINGLE_
#define _SE_GRID_AVX_256_SINGLE_

/* SE AVX 256 single gridding */
#if GMX_DOUBLE==0
#if GMX_SIMD_X86_AVX_256

#include "se.h"

// -----------------------------------------------------------------------------
static void SE_grid_split_AVX(real* gmx_restrict grid, real* gmx_restrict q,
			      splinedata_t * gmx_restrict spline, 
			      const SE_FGG_params* gmx_restrict params,
			      const pme_atomcomm_t * gmx_restrict atc,
			      const pmegrid_t * gmx_restrict pmegrid)
{
  // unpack parameters
  const int     N = params->N;
  float*        H = (float*) grid; // pointer to grid does NOT alias
  const float* zs = (float*) spline->zs;
  const float* zx = (float*) spline->theta[0];
  const float* zy = (float*) spline->theta[1];
  const float* zz = (float*) spline->theta[2];
    
  const int p = 8;
  float qn;
  int idx0, idx_zs, idx_zz, i, j, n, nn;
  int index_x, index_xy;
  int i0,j0,k0;
  int * idxptr;

  int pny = pmegrid->s[YY];
  int pnz = pmegrid->s[ZZ];

  int offx = pmegrid->offset[XX];
  int offy = pmegrid->offset[YY];
  int offz = pmegrid->offset[ZZ];

  __m256 rH0, rZZ0, rZS0, rC;

  for(n=0; n<N; n++){
    nn = spline->ind[n];
    qn = q[nn];
    idxptr = atc->idx[nn];
    i0 = idxptr[XX] - offx;
    j0 = idxptr[YY] - offy;
    k0 = idxptr[ZZ] - offz;
    idx0 = i0*pny*pnz+j0*pnz+k0;

    idx_zs = 0;

    //    if(idx0%8 == 0){ // H[idx0] is 32-aligned
    if(0){
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz;
	  rC = _mm256_set1_ps( qn*zx[p*n+i]*zy[p*n+j] );
	  idx_zz=p*n;
	  rH0  = _mm256_load_ps( H+index_xy + k0);
	  rZZ0 = _mm256_load_ps( zz + idx_zz     );
	  rZS0 = _mm256_load_ps( zs + idx_zs    );
	  
	  rZZ0 = _mm256_mul_ps(rZZ0,rC);
	  rZZ0 = _mm256_mul_ps(rZZ0,rZS0);
	  rH0  = _mm256_add_ps(rH0,rZZ0);
	  
	  _mm256_store_ps( H+index_xy + k0, rH0 );
	  
	  idx_zs+=8; 
	  idx_zz+=8;
	}
      }
    }
    else{ // H[idx0] is 16-aligned, preventing nice vectorization
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz;
	  rC = _mm256_set1_ps( qn*zx[p*n+i]*zy[p*n+j] );
	  idx_zz=p*n;
	  rH0  = _mm256_loadu_ps( H+index_xy + k0);

	  rZZ0 = _mm256_load_ps( zz + idx_zz );
	  rZS0 = _mm256_load_ps( zs + idx_zs );
	  
	  rZZ0 = _mm256_mul_ps(rZZ0,rC);
	  rZZ0 = _mm256_mul_ps(rZZ0,rZS0);
	  
	  rH0  = _mm256_add_ps(rH0,rZZ0);
	  _mm256_storeu_ps( H+index_xy + k0, rH0 );

	  idx_zs+=8;
	  idx_zz+=8;
	}
      }
    }
  }
}

// -----------------------------------------------------------------------------
static void 
SE_grid_split_AVX_dispatch(real* grid, real* q, 
			   splinedata_t *spline,
			   const SE_FGG_params* params,
			   const pme_atomcomm_t *atc,
			   const pmegrid_t *pmegrid)
{
  const int p = params->P;
  const int incrj = params->dims[2]; // middle increment
  const int incri = params->npdims[2]*(params->dims[1]);// outer increment

  // if P, or either increments are not divisible by 4, fall back on vanilla
  if( isnot_div_by_8(p) || isnot_div_by_8(incri) || isnot_div_by_8(incrj) || (p%8)!=0)
    {
      __DISPATCHER_MSG("[FGG GRID AVX SINGLE] AVX Abort (PARAMS)\n");
      SE_grid_split_SSE_dispatch(grid, q, spline, params,atc,pmegrid);
      return;
    }
    
  // otherwise the preconditions for AVX codes are satisfied. 
  else if(p==8){
    // specific for p=8
    __DISPATCHER_MSG("[FGG GRID AVX SINGLE] P=8\n");
    SE_grid_split_AVX(grid, q, spline, params, atc, pmegrid);
  }
}

#endif //GMX_SIMD_X86_AVX_256
#endif //not GMX_DOUBLE
#endif // _SE_GRID_AVX_256_SINGLE_
