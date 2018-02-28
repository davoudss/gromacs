#ifndef _SE_GRID_SSE_DOUBLE_
#define _SE_GRID_SSE_DOUBLE_

/* SE grid SSE double gridding */
#if GMX_DOUBLE==1
#include "se.h"

static void SE_grid_split_d(real* gmx_restrict grid, real* gmx_restrict q,
			      splinedata_t         * gmx_restrict spline, 
			      const pme_atomcomm_t * gmx_restrict atc,
			      const pmegrid_t      * gmx_restrict pmegrid) 
{
  // unpack parameters
  const real*   zs = (real*) spline->zs;
  const real*   zx = (real*) spline->theta[0];
  const real*   zy = (real*) spline->theta[1];
  const real*   zz = (real*) spline->theta[2];

  const int p = pmegrid->order-1;

  real cij0,qn;
  int idx0, zidx, idxzz, i, j, k, n, nn;
  int index_x, index_xy;
  int i0,j0,k0;
  int * idxptr;

  int pny = pmegrid->s[YY];
  int pnz = pmegrid->s[ZZ];

  int offx = pmegrid->offset[XX];
  int offy = pmegrid->offset[YY];
  int offz = pmegrid->offset[ZZ];

/*
#ifdef _OPENMP
#pragma omp for private(n) schedule(static) // work-share over OpenMP threads here
#endif
*/
  for(n=0; n<spline->n; n++)
    {
      nn = spline->ind[n];
      qn = q[nn];
      idxptr = atc->idx[nn];
      i0 = idxptr[XX] - offx;
      j0 = idxptr[YY] - offy;
      k0 = idxptr[ZZ] - offz;

      // inline vanilla loop
      zidx = 0;
      for(i = 0; i<p; i++)
        {
	  index_x = (i0+i)*pny*pnz;
	  real qnzx = qn*zx[p*n+i];
          for(j = 0; j<p; j++)
            {
	      index_xy = index_x + (j0+j)*pnz + k0;
              cij0 = qnzx*zy[p*n+j];
              idxzz=p*n;
              for(k = 0; k<p; k++)
                {
		  idx0 = index_xy + k;
                  grid[idx0] += zs[zidx]*zz[idxzz]*cij0;
                  zidx++; idxzz++;
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
static void SE_grid_split_SSE_d(real* gmx_restrict grid, real* gmx_restrict q,
				splinedata_t         * gmx_restrict spline, 
				const pme_atomcomm_t * gmx_restrict atc,
				const pmegrid_t      * gmx_restrict pmegrid)
{
  // unpack parameters
  const double*   zs = (double*) spline->zs;
  const double*   zx = (double*) spline->theta[0];
  const double*   zy = (double*) spline->theta[1];
  const double*   zz = (double*) spline->theta[2];
    
  const int p = pmegrid->order-1;
  double qn;
  int idx_zs, idx_zz, i, j, k, n, nn;
  int index_x, index_xy;
  int i0,j0,k0;
  int * idxptr;

  int pny = pmegrid->s[YY];
  int pnz = pmegrid->s[ZZ];

  int offx = pmegrid->offset[XX];
  int offy = pmegrid->offset[YY];
  int offz = pmegrid->offset[ZZ];

  __m128d rH0, rZZ0, rZS0, rC;

  for(n=0; n<spline->n; n++){
    nn = spline->ind[n];
    qn = q[nn];
    idxptr = atc->idx[nn];
    i0 = idxptr[XX] - offx;
    j0 = idxptr[YY] - offy;
    k0 = idxptr[ZZ] - offz;

    idx_zs = 0;
    
    if(0){ // H[idx0] is 16-aligned
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	real qnzx = qn*zx[p*n+i];
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  rC = _mm_set1_pd( qnzx * zy[p*n+j] );
	  idx_zz=p*n;
	  for(k = 0; k<p; k+=2){
	    rH0  = _mm_load_pd( grid + index_xy + k );
	    rZZ0 = _mm_load_pd( zz + idx_zz     );
	    rZS0 = _mm_load_pd( zs + idx_zs    );
	    rZZ0 = _mm_mul_pd(rZZ0,rC);
	    rZZ0 = _mm_mul_pd(rZZ0,rZS0);
	    rH0  = _mm_add_pd(rH0,rZZ0);
	    _mm_store_pd( grid + index_xy + k , rH0 );
	    
	    idx_zs+=2; 
	    idx_zz+=2;
	  }
	}
      }
    }
    else{ // H[idx0] is 8-aligned, preventing nice vectorization
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	real qnzx = qn*zx[p*n+i];
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  rC = _mm_set1_pd( qnzx * zy[p*n+j] );
	  idx_zz=p*n;
	  for(k = 0; k<p; k+=2){
	    rH0  = _mm_loadu_pd( grid + index_xy + k );
	    rZZ0 = _mm_load_pd( zz + idx_zz     );
	    rZS0 = _mm_load_pd( zs + idx_zs    );
	    rZZ0 = _mm_mul_pd(rZZ0,rC);
	    rZZ0 = _mm_mul_pd(rZZ0,rZS0);
	    rH0  = _mm_add_pd(rH0,rZZ0);
	    _mm_storeu_pd( grid + index_xy + k , rH0 );
	    
	    idx_zs+=2; 
	    idx_zz+=2;
	  }
	}
      }      
    }
  }
}

// -----------------------------------------------------------------------------
static void SE_grid_split_SSE_u8_d(real* gmx_restrict grid, real* gmx_restrict q,
				   splinedata_t         * gmx_restrict spline, 
				   const pme_atomcomm_t * gmx_restrict atc,
				   const pmegrid_t      * gmx_restrict pmegrid)
{
  // unpack parameters
  const double*   zs = (double*) spline->zs;
  const double*   zx = (double*) spline->theta[0];
  const double*   zy = (double*) spline->theta[1];
  const double*   zz = (double*) spline->theta[2];
    
  const int p = pmegrid->order-1;

  double qn;
  int idx0, idx_zs, idx_zz, i, j, k, n, nn;
  int index_x, index_xy;
  int i0,j0,k0;
  int * idxptr;

  int pny = pmegrid->s[YY];
  int pnz = pmegrid->s[ZZ];

  int offx = pmegrid->offset[XX];
  int offy = pmegrid->offset[YY];
  int offz = pmegrid->offset[ZZ];

  __m128d rH0, rZZ0, rZS0, rC;
  __m128d rH1, rZZ1, rZS1;
  __m128d rH2, rZZ2, rZS2;
  __m128d rH3, rZZ3, rZS3;

  for(n=0; n<spline->n; n++){
    nn = spline->ind[n];
    qn = q[nn];
    idxptr = atc->idx[nn];
    i0 = idxptr[XX] - offx;
    j0 = idxptr[YY] - offy;
    k0 = idxptr[ZZ] - offz;
    idx0 = i0*pny*pnz+j0*pnz+k0;
    _mm_prefetch( (void*) (grid + idx0), _MM_HINT_T0);
    
    idx_zs = 0;
    _mm_prefetch( (void*) zs, _MM_HINT_T0);
    
    if(0){ // H[idx0] is 16-aligned
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	real qnzx = qn*zx[p*n+i];
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  rC = _mm_set1_pd( qnzx*zy[p*n+j] );
	  idx_zz=p*n;
	  
	  for(k = 0; k<p; k+=8){
	    rH0  = _mm_load_pd( grid + index_xy + k );
	    rH1  = _mm_load_pd( grid + index_xy + k + 2 );
	    rH2  = _mm_load_pd( grid + index_xy + k + 4 );
	    rH3  = _mm_load_pd( grid + index_xy + k + 6 );
	    
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
	    
	    _mm_store_pd( grid + index_xy + k    , rH0 );
	    _mm_store_pd( grid + index_xy + k + 2, rH1 );
	    _mm_store_pd( grid + index_xy + k + 4, rH2 );
	    _mm_store_pd( grid + index_xy + k + 6, rH3 );
	    
	    idx_zs+=8; 
	    idx_zz+=8;
	  }
	}
      }
    }
    else{ // H[idx0] is 8-aligned, preventing nice vectorization
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	real qnzx = qn*zx[p*n+i];
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  rC = _mm_set1_pd( qnzx*zy[p*n+j] );
	  idx_zz=p*n;
	  
	  for(k = 0; k<p; k+=8){
	    rH0  = _mm_loadu_pd( grid + index_xy + k );
	    rH1  = _mm_loadu_pd( grid + index_xy + k + 2 );
	    rH2  = _mm_loadu_pd( grid + index_xy + k + 4 );
	    rH3  = _mm_loadu_pd( grid + index_xy + k + 6 );
	    
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
	    
	    _mm_storeu_pd( grid + index_xy + k    , rH0 );
	    _mm_storeu_pd( grid + index_xy + k + 2, rH1 );
	    _mm_storeu_pd( grid + index_xy + k + 4, rH2 );
	    _mm_storeu_pd( grid + index_xy + k + 6, rH3 );
	    
	    idx_zs+=8; 
	    idx_zz+=8;
	  }
	}
      }
    }
  }
}

// -----------------------------------------------------------------------------
static void SE_grid_split_SSE_P8_d(real* gmx_restrict grid, real* gmx_restrict q,
				   splinedata_t         * gmx_restrict spline, 
				   const pme_atomcomm_t * gmx_restrict atc,
				   const pmegrid_t      * gmx_restrict pmegrid)
{
  // unpack parameters
  const double*   zs = (double*) spline->zs;
  const double*   zx = (double*) spline->theta[0];
  const double*   zy = (double*) spline->theta[1];
  const double*   zz = (double*) spline->theta[2];
    
  double qn;
  int idx0, idx_zs, idx_zz, i, j, k, n, nn;
  int index_x, index_xy;
  int i0,j0,k0;
  int * idxptr;

  int pny = pmegrid->s[YY];
  int pnz = pmegrid->s[ZZ];

  int offx = pmegrid->offset[XX];
  int offy = pmegrid->offset[YY];
  int offz = pmegrid->offset[ZZ];

  __m128d rH0, rZZ0, rZS0, rC;
  __m128d rH1, rZZ1, rZS1;
  __m128d rH2, rZZ2, rZS2;
  __m128d rH3, rZZ3, rZS3;

  for(n=0; n<spline->n; n++){
    nn = spline->ind[n];
    qn = q[nn];
    idxptr = atc->idx[nn];
    i0 = idxptr[XX] - offx;
    j0 = idxptr[YY] - offy;
    k0 = idxptr[ZZ] - offz;
    idx0 = i0*pny*pnz+j0*pnz+k0;
    _mm_prefetch( (void*) (grid+idx0), _MM_HINT_T0);
    
    idx_zs = 0;
    _mm_prefetch( (void*) zs, _MM_HINT_T0);
    
    if(0){ // H[idx0] is 16-aligned
      for(i = 0; i<8; i++){
	index_x = (i0+i)*pny*pnz;
	real qnzx = qn*zx[8*n+i];
	for(j = 0; j<8; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  rC = _mm_set1_pd( qnzx * zy[8*n+j] );
	  idx_zz=8*n;
	  
	  for(k = 0; k<8; k+=8){
	    rH0  = _mm_load_pd( grid + index_xy + k );
	    rH1  = _mm_load_pd( grid + index_xy + k + 2 );
	    rH2  = _mm_load_pd( grid + index_xy + k + 4 );
	    rH3  = _mm_load_pd( grid + index_xy + k + 6 );

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
			
	    _mm_store_pd( grid + index_xy + k    , rH0 );
	    _mm_store_pd( grid + index_xy + k + 2, rH1 );
	    _mm_store_pd( grid + index_xy + k + 4, rH2 );
	    _mm_store_pd( grid + index_xy + k + 6, rH3 );

	    idx_zs+=8; 
	    idx_zz+=8;
	  }
	}
      }
    }
    else{ // H[idx0] is 8-aligned, preventing nice vectorization
      for(i = 0; i<8; i++){
	index_x = (i0+i)*pny*pnz;
	real qnzx = qn*zx[8*n+i];
	for(j = 0; j<8; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  rC = _mm_set1_pd( qnzx * zy[8*n+j] );
	  idx_zz=8*n;
	  
	  for(k = 0; k<8; k+=8){
	    rH0  = _mm_loadu_pd( grid + index_xy + k );
	    rH1  = _mm_loadu_pd( grid + index_xy + k + 2 );
	    rH2  = _mm_loadu_pd( grid + index_xy + k + 4 );
	    rH3  = _mm_loadu_pd( grid + index_xy + k + 6 );

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
			
	    _mm_storeu_pd( grid + index_xy + k    , rH0 );
	    _mm_storeu_pd( grid + index_xy + k + 2, rH1 );
	    _mm_storeu_pd( grid + index_xy + k + 4, rH2 );
	    _mm_storeu_pd( grid + index_xy + k + 6, rH3 );

	    idx_zs+=8; 
	    idx_zz+=8;
	  }
	}
      }
    }
  }
}


// -----------------------------------------------------------------------------
static void SE_grid_split_SSE_P16_d(real* gmx_restrict grid, real* gmx_restrict q,
				    splinedata_t         * gmx_restrict spline, 
				    const pme_atomcomm_t * gmx_restrict atc,
				    const pmegrid_t      * gmx_restrict pmegrid)
{
  // unpack parameters
  const double* zs = (double*) spline->zs;
  const double* zx = (double*) spline->theta[0];
  const double* zy = (double*) spline->theta[1];
  const double* zz = (double*) spline->theta[2];

  int idx0, idx_zs, i, j, n, nn;
  double qn;
  int index_x, index_xy;
  int i0,j0,k0;
  int * idxptr;

  int pny = pmegrid->s[YY];
  int pnz = pmegrid->s[ZZ];

  int offx = pmegrid->offset[XX];
  int offy = pmegrid->offset[YY];
  int offz = pmegrid->offset[ZZ];


  __m128d rZZ0, rZZ1, rZZ2, rZZ3, rZZ4, rZZ5, rZZ6, rZZ7; 
  __m128d rH0, rH1, rH2, rH3;
  __m128d rC, rZS0;

  for(n=0; n<spline->n; n++){
    nn = spline->ind[n];
    qn = q[nn];
    idxptr = atc->idx[nn];
    i0 = idxptr[XX] - offx;
    j0 = idxptr[YY] - offy;
    k0 = idxptr[ZZ] - offz;
    idx0 = i0*pny*pnz+j0*pnz+k0;
    _mm_prefetch( (void*) (grid+idx0), _MM_HINT_T0);

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

    if(0){ // H[idx0] is 16-aligned
      for(i = 0; i<16; i++){
	index_x = (i0+i)*pny*pnz;
	real qnzx = qn*zx[16*n+i];
	for(j = 0; j<16; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  rC = _mm_set1_pd( qnzx * zy[16*n+j] );

	  /* 0 - 3 */ 
	  rH0  = _mm_load_pd( grid + index_xy    );
	  rH1  = _mm_load_pd( grid + index_xy + 2);
	  rH2  = _mm_load_pd( grid + index_xy + 4);
	  rH3  = _mm_load_pd( grid + index_xy + 6);

	  rZS0 = _mm_load_pd( zs + idx_zs);
	  rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0));

	  rZS0 = _mm_load_pd( zs + idx_zs + 2);                   
	  rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS0));

	  rZS0 = _mm_load_pd( zs + idx_zs + 4);                   
	  rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS0));

	  rZS0 = _mm_load_pd( zs + idx_zs + 6);                   
	  rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS0));

	  _mm_store_pd(grid + index_xy    , rH0);
	  _mm_store_pd(grid + index_xy + 2, rH1);
	  _mm_store_pd(grid + index_xy + 4, rH2);
	  _mm_store_pd(grid + index_xy + 6, rH3);

	  /* 4 - 7*/ 
	  rH0  = _mm_load_pd( grid + index_xy + 8 );
	  rH1  = _mm_load_pd( grid + index_xy + 10);
	  rH2  = _mm_load_pd( grid + index_xy + 12);
	  rH3  = _mm_load_pd( grid + index_xy + 14);

	  rZS0 = _mm_load_pd( zs + idx_zs + 8);
	  rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rC),rZS0));

	  rZS0 = _mm_load_pd( zs + idx_zs + 10);                   
	  rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ5,rC),rZS0));

	  rZS0 = _mm_load_pd( zs + idx_zs + 12);                   
	  rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ6,rC),rZS0));

	  rZS0 = _mm_load_pd( zs + idx_zs + 14);                   
	  rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ7,rC),rZS0));

	  _mm_store_pd(grid + index_xy + 8 , rH0);
	  _mm_store_pd(grid + index_xy + 10, rH1);
	  _mm_store_pd(grid + index_xy + 12, rH2);
	  _mm_store_pd(grid + index_xy + 14, rH3);

	  idx_zs += 16;
	}
      }
    }
    else{ // H[idx0] is 8-aligned, preventing nice vectorization
      for(i = 0; i<16; i++){
	index_x = (i0+i)*pny*pnz;
	real qnzx = qn*zx[16*n+i];
	for(j = 0; j<16; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  rC = _mm_set1_pd( qnzx * zy[16*n+j] );

	  /* 0 - 3 */ 
	  rH0  = _mm_loadu_pd( grid + index_xy    );
	  rH1  = _mm_loadu_pd( grid + index_xy + 2);
	  rH2  = _mm_loadu_pd( grid + index_xy + 4);
	  rH3  = _mm_loadu_pd( grid + index_xy + 6);

	  rZS0 = _mm_load_pd( zs + idx_zs);
	  rH0  = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0));

	  rZS0 = _mm_load_pd( zs + idx_zs + 2);                   
	  rH1  = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS0));

	  rZS0 = _mm_load_pd( zs + idx_zs + 4);                   
	  rH2  = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS0));

	  rZS0 = _mm_load_pd( zs + idx_zs + 6);                   
	  rH3  = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS0));

	  _mm_storeu_pd(grid + index_xy    , rH0);
	  _mm_storeu_pd(grid + index_xy + 2, rH1);
	  _mm_storeu_pd(grid + index_xy + 4, rH2);
	  _mm_storeu_pd(grid + index_xy + 6, rH3);

	  /* 4 - 7*/ 
	  rH0  = _mm_loadu_pd( grid + index_xy + 8 );
	  rH1  = _mm_loadu_pd( grid + index_xy + 10);
	  rH2  = _mm_loadu_pd( grid + index_xy + 12);
	  rH3  = _mm_loadu_pd( grid + index_xy + 14);

	  rZS0 = _mm_load_pd( zs + idx_zs + 8);
	  rH0  = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rC),rZS0));

	  rZS0 = _mm_load_pd( zs + idx_zs + 10);                   
	  rH1  = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ5,rC),rZS0));

	  rZS0 = _mm_load_pd( zs + idx_zs + 12);                   
	  rH2  = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ6,rC),rZS0));

	  rZS0 = _mm_load_pd( zs + idx_zs + 14);                   
	  rH3  = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ7,rC),rZS0));

	  _mm_storeu_pd(grid + index_xy + 8 , rH0);
	  _mm_storeu_pd(grid + index_xy + 10, rH1);
	  _mm_storeu_pd(grid + index_xy + 12, rH2);
	  _mm_storeu_pd(grid + index_xy + 14, rH3);

	  idx_zs += 16;
	}
      }      
    }
  }
}

// -----------------------------------------------------------------------------
static void
SE_grid_split_SSE_dispatch_d(real* grid, real* q, 
			   splinedata_t *spline,
			   const SE_params* params,
			   const pme_atomcomm_t *atc,
			   const pmegrid_t *pmegrid)
{
  const int p = params->P;
  const int incrj = params->dims[2]; // middle increment
  const int incri = params->npdims[2]*(params->dims[1]);// outer increment

  // if P is odd, or if either increment is odd, fall back on vanilla
  if( is_odd(p) || is_odd(incri) || is_odd(incrj) ){
    __DISPATCHER_MSG("[FGG GRID SSE DOUBLE] SSE Abort (PARAMS)\n");
    SE_grid_split_d(grid, q, spline, atc, pmegrid);
    return;
  }

  // otherwise the preconditions for SSE codes are satisfied.  
  if(p==16){
    // specific for p=16
    __DISPATCHER_MSG("[FGG GRID SSE DOUBLE] P=16\n");
    SE_grid_split_SSE_P16_d(grid, q, spline,atc, pmegrid);
  }
  else if(p==8){
    // specific for p divisible by 8
    __DISPATCHER_MSG("[FGG GRID SSE DOUBLE] P=8\n");
    SE_grid_split_SSE_P8_d(grid, q, spline, atc, pmegrid); 
  }  
  else if(p%8==0){
    // specific for p divisible by 8
    __DISPATCHER_MSG("[FGG GRID SSE DOUBLE] P unroll 8\n");
    SE_grid_split_SSE_u8_d(grid, q, spline, atc, pmegrid); 
  }
  else{
    // vanilla SSE code (any even p)
    __DISPATCHER_MSG("[FGG GRID SSE DOUBLE] Vanilla\n");
    SE_grid_split_SSE_d(grid, q, spline, atc, pmegrid);
  }
}

#endif // GMX_DOUBLE
#endif //_SE_GRID_SSE_DOUBLE_
