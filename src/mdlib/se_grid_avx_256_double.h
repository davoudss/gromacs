#ifndef _SE_GRID_AVX_256_DOUBLE_
#define _SE_GRID_AVX_256_DOUBLE_


/* SE AVX 256 double gridding */
#ifdef GMX_DOUBLE
#ifdef GMX_X86_AVX_256
#include "se.h"


// -----------------------------------------------------------------------------
void SE_grid_split_AVX_d(real* grid, real* q,
			 splinedata_t *spline, 
			 const SE_FGG_params* params)
{
  // unpack parameters
  const int      N = params->N;
  double*        H = (double*) grid; // pointer to grid does NOT alias
  const double* zs = (double*) spline->zs;
  const double* zx = (double*) spline->theta[0];
  const double* zy = (double*) spline->theta[1];
  const double* zz = (double*) spline->theta[2];
    
  const int p = params->P;
  double qn;
  int idx0, idx_zs, idx_zz, i, j, k, n, nn;
  const int incrj = params->npdims[2]-p; // middle increment
  const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

  __m256d rH0, rZZ0, rZS0, rC;

  for(n=0; n<N; n++)
    {
      nn = spline->ind[n];
      qn = q[nn];
      idx0 = spline->idx[n];
      idx_zs = 0;
      
      if(idx0%4 == 0) // H[idx0] is 32-aligned
	{
	  for(i = 0; i<p; i++)
	    {
	      for(j = 0; j<p; j++)
		{
		  rC = _mm256_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
		  idx_zz=p*n;
		  for(k = 0; k<p; k+=4)
		    {
		      rH0  = _mm256_load_pd( H+idx0     );
		      rZZ0 = _mm256_load_pd( zz + idx_zz     );
		      rZS0 = _mm256_load_pd( zs + idx_zs    );

		      rZZ0 = _mm256_mul_pd(rZZ0,rC);
		      rZZ0 = _mm256_mul_pd(rZZ0,rZS0);
		      rH0  = _mm256_add_pd(rH0,rZZ0);

		      _mm256_store_pd( H+idx0    , rH0 );

		      idx0  +=4;
		      idx_zs+=4; 
		      idx_zz+=4;
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
		  rC = _mm256_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
		  idx_zz=p*n;
		  for(k = 0; k<p; k+=4)
	    	    {
		      rH0  = _mm256_loadu_pd( H+idx0 );

		      rZZ0 = _mm256_load_pd( zz + idx_zz );
		      rZS0 = _mm256_load_pd( zs + idx_zs );

		      rZZ0 = _mm256_mul_pd(rZZ0,rC);
		      rZZ0 = _mm256_mul_pd(rZZ0,rZS0);

		      rH0  = _mm256_add_pd(rH0,rZZ0);
		      _mm256_storeu_pd( H+idx0, rH0 );

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
void SE_grid_split_AVX_u8_d(real* grid, real* q,
			    splinedata_t *spline,
			    const SE_FGG_params* params)
{
  // unpack parameters
  const int      N = params->N;
  double*        H = (double*) grid; // pointer to grid does NOT alias
  const double* zs = (double*) spline->zs;
  const double* zx = (double*) spline->theta[0];
  const double* zy = (double*) spline->theta[1];
  const double* zz = (double*) spline->theta[2];
    
  const int p = params->P;

  const int incrj = params->npdims[2]-p; // middle increment
  const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

  double qn;
  int idx0, idx_zs, idx_zz, i, j, k, n, nn;

  __m256d rH0, rZZ0, rZS0, rC;
  __m256d rH1, rZZ1, rZS1;

  for(n=0; n<N; n++){
    nn = spline->ind[n];
    qn = q[nn];
    idx0 = spline->idx[n];
    _mm_prefetch( (void*) (H+idx0), _MM_HINT_T0); 
    idx_zs = 0;
    _mm_prefetch( (void*) zs, _MM_HINT_T0);
    
    if(idx0%4 == 0){ // H[idx0] is 32-aligned
      for(i = 0; i<p; i++){
	for(j = 0; j<p; j++){
	  rC = _mm256_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
	  idx_zz=p*n;
	  
	  for(k = 0; k<p; k+=8){
	    rH0  = _mm256_load_pd( H+idx0     );
	    rH1  = _mm256_load_pd( H+idx0 + 4 );
		      
	    rZZ0 = _mm256_load_pd( zz + idx_zz     );
	    rZZ1 = _mm256_load_pd( zz + idx_zz + 4 );
		      
	    rZS0 = _mm256_load_pd( zs + idx_zs    );
	    rZS1 = _mm256_load_pd( zs + idx_zs + 4);
		      
	    rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
	    rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
		      
	    _mm256_store_pd( H+idx0    , rH0 );
	    _mm256_store_pd( H+idx0 + 4, rH1 );
		      
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
	  rC = _mm256_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
	  idx_zz=p*n;
	  
	  for(k = 0; k<p; k+=8){
	    rH0  = _mm256_loadu_pd( H+idx0     );
	    rH1  = _mm256_loadu_pd( H+idx0 + 4 );
		      
	    rZZ0 = _mm256_load_pd( zz + idx_zz     );
	    rZZ1 = _mm256_load_pd( zz + idx_zz + 4 );
		      
	    rZS0 = _mm256_load_pd( zs + idx_zs    );
	    rZS1 = _mm256_load_pd( zs + idx_zs + 4);
		      
	    rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
	    rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
		      
	    _mm256_storeu_pd( H+idx0    , rH0 );
	    _mm256_storeu_pd( H+idx0 + 4, rH1 );
		      
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
void SE_grid_split_AVX_P8_d(real* grid, real* q,
			    splinedata_t *spline,
			    const SE_FGG_params* params)
{
  // unpack parameters
  const int      N = params->N;
  double*        H = (double*) grid; // pointer to grid does NOT alias
  const double* zs = (double*) spline->zs;
  const double* zx = (double*) spline->theta[0];
  const double* zy = (double*) spline->theta[1];
  const double* zz = (double*) spline->theta[2];

  double qn;
  int idx, idx_zs, i, j, n, nn;
  const int incrj = params->npdims[2]-8; // middle increment
  const int incri = params->npdims[2]*(params->npdims[1]-8);// outer increment
   
  __m256d rZZ0, rZZ1; 
  __m256d rH0, rH1;
  __m256d rC, rZS0,rZS1;

  for(n=0; n<spline->n; n++){
    nn = spline->ind[n];
    qn = q[nn];
    idx = spline->idx[n];
    idx_zs = 0;
    
    rZZ0 = _mm256_load_pd(zz + n*8     );
    rZZ1 = _mm256_load_pd(zz + n*8 + 4 );
    
    if(idx%4 == 0){ // H[idx0] is 32-aligned
      for(i = 0; i<8; i++){
	for(j = 0; j<8; j++){
	  rC = _mm256_set1_pd( qn*zx[8*n+i]*zy[8*n+j] );

	  rH0  = _mm256_load_pd( H+idx    );
	  rH1  = _mm256_load_pd( H+idx + 4);

	  rZS0 = _mm256_load_pd( zs + idx_zs);
	  rZS1 = _mm256_load_pd( zs + idx_zs + 4);

	  rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
	  rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
		  
	  _mm256_store_pd(H + idx,      rH0);
	  _mm256_store_pd(H + idx + 4,  rH1);

	  idx += incrj + 8;
	  idx_zs += 8;
	}
	idx += incri;
      }
    }
    else{ // H[idx0] is 16-aligned, preventing nice vectorization
      for(i = 0; i<8; i++){
	for(j = 0; j<8; j++){
	  rC = _mm256_set1_pd( qn*zx[8*n+i]*zy[8*n+j] );
		  
	  rH0  = _mm256_loadu_pd( H+idx     );
	  rH1  = _mm256_loadu_pd( H+idx + 4 );
		  
	  rZS0 = _mm256_load_pd( zs + idx_zs     );
	  rZS1 = _mm256_load_pd( zs + idx_zs + 4 );
		  
	  rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
	  rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
		  
	  _mm256_storeu_pd(H + idx,      rH0);
	  _mm256_storeu_pd(H + idx + 4,  rH1);
		  
	  idx += incrj + 8;
	  idx_zs += 8;
	}
	idx += incri;
      }
    }
  }
}

// -------------------------------------------------------------------------
void SE_grid_split_AVX_P16_d(real *grid, real* q,
			     splinedata_t *spline,
			     const SE_FGG_params* params)
{
  // unpack parameters
  const int      N = params->N;
  double*        H = (double*) grid; // pointer to grid does NOT alias
  const double* zs = (double*) spline->zs;
  const double* zx = (double*) spline->theta[0];
  const double* zy = (double*) spline->theta[1];
  const double* zz = (double*) spline->theta[2];


  const int incrj = params->npdims[2]-16; // middle increment
  const int incri = params->npdims[2]*(params->npdims[1]-16);// outer increment

  double qn;
  int idx, idx_zs, i, j, n, nn;

  __m256d rZZ0, rZZ1, rZZ2, rZZ3; 
  __m256d rH0, rH1, rH2, rH3;
  __m256d rC, rZS0,rZS1,rZS2,rZS3;

  for(n=0; n<N; n++){
    nn = spline->ind[n];
    qn = q[nn];
    idx = spline->idx[n];
    _mm_prefetch( (void*) (H+idx), _MM_HINT_T0);
    idx_zs = 0;
    _mm_prefetch( (void*) zs, _MM_HINT_T0);

    rZZ0 = _mm256_load_pd(zz + n*16     );
    rZZ1 = _mm256_load_pd(zz + n*16 + 4 );
    rZZ2 = _mm256_load_pd(zz + n*16 + 8 );
    rZZ3 = _mm256_load_pd(zz + n*16 + 12);
      
    if(idx%4 == 0){ // H[idx0] is 32-aligned
      for(i = 0; i<16; i++){
	for(j = 0; j<16; j++){
	  rC = _mm256_set1_pd( qn*zx[16*n+i]*zy[16*n+j] );
		  
	  rH0  = _mm256_load_pd( H+idx    );
	  rH1  = _mm256_load_pd( H+idx + 4);
	  rH2  = _mm256_load_pd( H+idx + 8);
	  rH3  = _mm256_load_pd( H+idx + 12);

	  rZS0 = _mm256_load_pd( zs + idx_zs);
	  rZS1 = _mm256_load_pd( zs + idx_zs + 4);
	  rZS2 = _mm256_load_pd( zs + idx_zs + 8);   
	  rZS3 = _mm256_load_pd( zs + idx_zs + 12);    

	  rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
	  rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
	  rH2 = _mm256_add_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rC),rZS2));
	  rH3 = _mm256_add_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rC),rZS3));

	  _mm256_store_pd(H + idx,      rH0);
	  _mm256_store_pd(H + idx + 4,  rH1);
	  _mm256_store_pd(H + idx + 8,  rH2);
	  _mm256_store_pd(H + idx + 12, rH3);

	  idx += incrj + 16;
	  idx_zs += 16;
	}
	idx += incri;
      }
    }
    else{ // H[idx0] is 8-aligned, preventing nice vectorization
      for(i = 0; i<16; i++){
	for(j = 0; j<16; j++){
	  rC = _mm256_set1_pd( qn*zx[16*n+i]*zy[16*n+j] );

	  rH0  = _mm256_loadu_pd( H+idx     );
	  rH1  = _mm256_loadu_pd( H+idx + 4 );
	  rH2  = _mm256_loadu_pd( H+idx + 8 );
	  rH3  = _mm256_loadu_pd( H+idx + 12);

	  rZS0 = _mm256_load_pd( zs + idx_zs     );
	  rZS1 = _mm256_load_pd( zs + idx_zs + 4 );
	  rZS2 = _mm256_load_pd( zs + idx_zs + 8 );
	  rZS3 = _mm256_load_pd( zs + idx_zs + 12);

	  rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
	  rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
	  rH2 = _mm256_add_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rC),rZS2));
	  rH3 = _mm256_add_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rC),rZS3));

	  _mm256_storeu_pd(H + idx,      rH0);
	  _mm256_storeu_pd(H + idx + 4,  rH1);
	  _mm256_storeu_pd(H + idx + 8,  rH2);
	  _mm256_storeu_pd(H + idx + 12, rH3);

	  idx += incrj + 16;
	  idx_zs += 16;
	}
	idx += incri;
      }
    }
  }
}

// ---------------------------------------------
void
SE_grid_split_AVX_dispatch_d(real* grid, real* q, 
			     splinedata_t *spline,
			     const SE_FGG_params* params)
{
  const int p = params->P;
  const int incrj = params->dims[2]; // middle increment
  const int incri = params->npdims[2]*(params->dims[1]);// outer increment
  
  // if P, or either increments are not divisible by 4,fall back on vanilla
  if( isnot_div_by_4(p) || isnot_div_by_4(incri) || isnot_div_by_4(incrj) ) {
    __DISPATCHER_MSG("[FGG GRID AVX] AVX Abort (PARAMS)\n");
    //	SE_grid_split_SSE_d(grid, q, spline, params);
    SE_grid_split_SSE_dispatch_d(grid, q, spline, params);
    return;
  }
  
  // otherwise the preconditions for AVX codes are satisfied. 
  if(p==16){
    // specific for p=16
    __DISPATCHER_MSG("[FGG GRID AVX] P=16\n");
    SE_grid_split_AVX_P16_d(grid, q, spline, params); 
  }
  else if(p==8){
    // specific for p=8
    __DISPATCHER_MSG("[FGG GRID AVX] P=8\n");
    SE_grid_split_AVX_P8_d(grid, q, spline, params); 
  }
  else if(p%8==0){
    // specific for p divisible by 8
    __DISPATCHER_MSG("[FGG GRID AVX] P unroll 8\n");
    SE_grid_split_AVX_u8_d(grid, q, spline, params); 
  }
  else if(p%4==0){
    // specific for p divisible by 4
    __DISPATCHER_MSG("[FGG GRID AVX] P unroll 4\n");
    SE_grid_split_AVX_d(grid, q, spline, params); 
  }
  else{
    // vanilla SSE code (any even p)
    __DISPATCHER_MSG("[FGG GRID AVX] Vanilla\n");
    SE_grid_split_SSE_d(grid, q, spline, params);
  } 
}


#endif //GMX_X86_AVX_256
#endif //GMX_DOUBLE
#endif //_SE_GRID_AVX_256_DOUBLE_
