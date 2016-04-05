#ifndef _SE_GRID_SSE_SINGLE_
#define _SE_GRID_SSE_SINGLE_

/* SE int SSE single gridding */
#ifndef GMX_DOUBLE 
#include "se.h"

// -----------------------------------------------------------------------------
void SE_grid_split(real* grid, real* q,
		   splinedata_t *spline,
		   const SE_FGG_params* params)
{
  // unpack parameters
  const int       N = params->N;
  float*          H = (float*) grid; // pointer to grid does NOT alias
  const float*   zs = (float*) spline->zs;
  const float*   zx = (float*) spline->theta[0];
  const float*   zy = (float*) spline->theta[1];
  const float*   zz = (float*) spline->theta[2];
    
  const int p = params->P;

  float cij0,qn;
  int idx0, zidx, idxzz, i, j, k, n, nn;
  const int incrj = params->npdims[2]-p; // middle increment
  const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

#ifdef _OPENMP
  // work-share over OpenMP threads here
#pragma omp for private(n) schedule(static)
#endif
  for(n=0; n<spline->n; n++)
    {
      nn = spline->ind[n];
      qn = q[nn];
      idx0 = spline->idx[n];

      // inline vanilla loop
      zidx = 0;
      for(i = 0; i<p; i++)
	{
	  for(j = 0; j<p; j++)
	    {
	      cij0 = zx[p*n+i]*zy[p*n+j];
	      idxzz=p*n;
	      for(k = 0; k<p; k++)
		{
		  H[idx0] += zs[zidx]*zz[idxzz]*cij0*qn;
		  idx0++; zidx++; idxzz++;
		}
	      idx0 += incrj; 
	    }
	  idx0 += incri; 
	}
    }
}


// -----------------------------------------------------------------------------
void SE_grid_split_SSE(real* grid, real* q,
		       splinedata_t *spline,
		       const SE_FGG_params* params)
{
  // unpack parameters
  const int       N = params->N;
  float*          H = (float*) grid; // pointer to grid does NOT alias
  const float*   zs = (float*) spline->zs;
  const float*   zx = (float*) spline->theta[0];
  const float*   zy = (float*) spline->theta[1];
  const float*   zz = (float*) spline->theta[2];
    
  const int p = params->P;

  int idx0, idx_zs, idx_zz, i, j, k, n, nn;
  float qn;
  const int incrj = params->npdims[2]-p; // middle increment
  const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

  __m128 rH0, rZZ0, rZS0, rC;

  for(n=0; n<spline->n; n++){
    nn = spline->ind[n];
    qn = q[nn];
    idx0 = spline->idx[n];
    idx_zs = 0;
    
    if(idx0%4 == 0){ // H[idx0] is 16-aligned
      for(i = 0; i<p; i++){
	for(j = 0; j<p; j++){
	  rC = _mm_set1_ps( qn*zx[p*n+i]*zy[p*n+j] );
	  idx_zz=p*n;
	  for(k = 0; k<p; k+=4){
	    rH0  = _mm_load_ps( H+idx0     );
	    rZZ0 = _mm_load_ps( zz + idx_zz     );
	    rZS0 = _mm_load_ps( zs + idx_zs    );
	    
	    rZZ0 = _mm_mul_ps(rZZ0,rC);
	    rZZ0 = _mm_mul_ps(rZZ0,rZS0);
	    rH0  = _mm_add_ps(rH0,rZZ0);
	    
	    _mm_store_ps( H+idx0    , rH0 );
	    
	    idx0  +=4;
	    idx_zs+=4; 
	    idx_zz+=4;
	  }
	  idx0 += incrj; 
	}
	idx0 += incri; 
      }
    }
    else{ // H[idx0] is 8-aligned, preventing nice vectorization
      for(i = 0; i<p; i++){
	for(j = 0; j<p; j++){
	  rC = _mm_set1_ps( qn*zx[p*n+i]*zy[p*n+j] );
	  idx_zz=p*n;
	  for(k = 0; k<p; k+=4){
	    rH0  = _mm_loadu_ps( H+idx0 );
	    rZZ0 = _mm_load_ps( zz + idx_zz );
	    rZS0 = _mm_load_ps( zs + idx_zs );
	    rZZ0 = _mm_mul_ps(rZZ0,rC);
	    rZZ0 = _mm_mul_ps(rZZ0,rZS0);
	    rH0  = _mm_add_ps(rH0,rZZ0);
	    _mm_storeu_ps( H+idx0, rH0 );
	    
	    idx0  +=4;
	    idx_zs+=4;
	    idx_zz+=4;
	  }
	  idx0 += incrj;
	}
	idx0 += incri;
      }
    }
  }
}

// -----------------------------------------------------------------------------
void SE_grid_split_SSE_P8(real *grid, real *q,
			  splinedata_t *spline,
			  const SE_FGG_params *params)
{
  // unpack parameters
  const int       N = params->N;
  float*          H = (float*) grid; // pointer to grid does NOT alias
  const float*   zs = (float*) spline->zs;
  const float*   zx = (float*) spline->theta[0];
  const float*   zy = (float*) spline->theta[1];
  const float*   zz = (float*) spline->theta[2];

  float qn;
  int idx, idx_zs, i, j, n , nn;
  const int incrj = params->npdims[2]-8; // middle increment
  const int incri = params->npdims[2]*(params->npdims[1]-8);// outer increment

  __m128 rZZ0, rZZ1; 
  __m128 rH0, rH1;
  __m128 rC, rZS0;

  for(n=0; n<spline->n; n++){
    nn = spline->ind[n];
    qn = q[nn];
    idx = spline->idx[n];
    _mm_prefetch( (void*) (H+idx), _MM_HINT_T0);
    
    idx_zs = 0;
    _mm_prefetch( (void*) zs, _MM_HINT_T0);
    
    rZZ0 = _mm_load_ps(zz + n*8     );
    rZZ1 = _mm_load_ps(zz + n*8 + 4 );
    if(idx%4 == 0){ // H[idx0] is 16-aligned
      for(i = 0; i<8; i++){
	for(j = 0; j<8; j++){
	  rC = _mm_set1_ps( qn*zx[8*n+i]*zy[8*n+j] );
	  
	  rH0  = _mm_load_ps( H+idx    );
	  rH1  = _mm_load_ps( H+idx + 4);
	  rZS0 = _mm_load_ps( zs + idx_zs);
	  rH0 = _mm_add_ps(rH0,_mm_mul_ps(_mm_mul_ps(rZZ0,rC),rZS0));
	  
	  rZS0 = _mm_load_ps( zs + idx_zs + 4);                   
	  rH1 = _mm_add_ps(rH1,_mm_mul_ps(_mm_mul_ps(rZZ1,rC),rZS0));
	  
	  _mm_store_ps(H + idx, rH0);
	  _mm_store_ps(H + idx + 4, rH1);
	  
	  idx += incrj + 8;
	  idx_zs += 8;
	}
	idx += incri;
      }
    }
    else{ // H[idx0] is 8-aligned, preventing nice vectorization
      for(i = 0; i<8; i++){
	for(j = 0; j<8; j++){
	  rC = _mm_set1_ps( qn*zx[8*n+i]*zy[8*n+j] );
	  
	  rH0  = _mm_loadu_ps( H+idx    );
	  rH1  = _mm_loadu_ps( H+idx + 4);
	  
	  // if zs does not have 16-byte alignment, this will core.
	  // PLATFORM AND COMPILER DEPENDENT (FIXME)
	  rZS0 = _mm_load_ps( zs + idx_zs);
	  rH0 = _mm_add_ps(rH0,_mm_mul_ps(_mm_mul_ps(rZZ0,rC),rZS0));
	  
	  rZS0 = _mm_load_ps( zs + idx_zs + 4);                   
	  rH1 = _mm_add_ps(rH1,_mm_mul_ps(_mm_mul_ps(rZZ1,rC),rZS0));
	  
	  _mm_storeu_ps(H + idx, rH0);
	  _mm_storeu_ps(H + idx + 4, rH1);
	  
	  idx += incrj + 8;
	  idx_zs += 8;
	}
	idx += incri;
      }
    }
  }
}

// -----------------------------------------------------------------------------
void 
SE_grid_split_SSE_dispatch(real* grid, real* q, 
			   splinedata_t *spline,
			   const SE_FGG_params* params)
{
  const int p = params->P;
  const int incrj = params->dims[2]; // middle increment
  const int incri = params->npdims[2]*(params->dims[1]);// outer increment

  // if P is odd, or if either increment is odd, fall back on vanilla
  if( isnot_div_by_4(p) || isnot_div_by_4(incri) || isnot_div_by_4(incrj)  || (p%4)!=0)
    {
      __DISPATCHER_MSG("[FGG GRID SSE] SSE Abort (PARAMS)\n");
      SE_grid_split(grid, q, spline, params);
      return;
    }
   
  // otherwise the preconditions for SSE codes are satisfied.    
  if(p==8){
    // specific for p=8
    __DISPATCHER_MSG("[FGG GRID SSE] P=8\n");
    SE_grid_split_SSE_P8(grid, q, spline, params);
  }
  else {
    // vanilla SSE code for p divisible by 4
    __DISPATCHER_MSG("[FGG GRID SSE] Vanilla\n");
    SE_grid_split_SSE(grid, q, spline, params);
  }
}

#endif // not GMX_DOUBLE

#endif //_SE_GRID_SSE_SINGLE_
