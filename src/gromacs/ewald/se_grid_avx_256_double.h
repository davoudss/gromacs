#ifndef _SE_GRID_AVX_256_DOUBLE_
#define _SE_GRID_AVX_256_DOUBLE_


/* SE AVX 256 double gridding */
#if GMX_DOUBLE==1
#if GMX_SIMD_X86_AVX_256
#include "se.h"


// -----------------------------------------------------------------------------
static
void SE_grid_split_AVX_gaussian_d(real* gmx_restrict grid, real* gmx_restrict q,
				  splinedata_t         * gmx_restrict spline, 
				  const pme_atomcomm_t * gmx_restrict atc,
				  const pmegrid_t      * gmx_restrict pmegrid)
{
  // unpack parameters
  const double* zs = (double*) spline->zs;
  const double* zx = (double*) spline->theta[0];
  const double* zy = (double*) spline->theta[1];
  const double* zz = (double*) spline->theta[2];
    
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

  __m256d rH0, rZZ0, rZS0, rC;

  for(n=0; n<spline->n; n++)
    {
      nn = spline->ind[n];
      qn = q[nn];
      idxptr = atc->idx[nn];
      i0 = idxptr[XX] - offx;
      j0 = idxptr[YY] - offy;
      k0 = idxptr[ZZ] - offz;
      idx_zs = 0;
      
      if(0) // H[idx0] is 32-aligned
	{
	  for(i = 0; i<p; i++)
	    {
	      index_x = (i0+i)*pny*pnz;
	      real qnzx = qn*zx[p*n+i];
	      for(j = 0; j<p; j++)
		{
		  index_xy = index_x + (j0+j)*pnz + k0;
		  rC = _mm256_set1_pd( qnzx*zy[p*n+j] );
		  idx_zz=p*n;
		  for(k = 0; k<p; k+=4)
		    {
		      rH0  = _mm256_load_pd( grid + index_xy + k );
		      rZZ0 = _mm256_load_pd( zz + idx_zz     );
		      rZS0 = _mm256_load_pd( zs + idx_zs    );

		      rZZ0 = _mm256_mul_pd(rZZ0,rC);
		      rZZ0 = _mm256_mul_pd(rZZ0,rZS0);
		      rH0  = _mm256_add_pd(rH0,rZZ0);

		      _mm256_store_pd(grid+index_xy + k , rH0 );

		      idx_zs+=4; 
		      idx_zz+=4;
		    }
		}
	    }
	}
      else // H[idx0] is 8-aligned, preventing nice vectorization
	{
	  for(i = 0; i<p; i++)
	    {
	      index_x = (i0+i)*pny*pnz;
	      real qnzx = qn*zx[p*n+i];
	      for(j = 0; j<p; j++)
		{
		  index_xy = index_x + (j0+j)*pnz + k0;
		  rC = _mm256_set1_pd( qnzx*zy[p*n+j] );
		  idx_zz=p*n;
		  for(k = 0; k<p; k+=4)
		    {
		      rH0  = _mm256_loadu_pd( grid + index_xy + k );
		      rZZ0 = _mm256_load_pd( zz + idx_zz     );
		      rZS0 = _mm256_load_pd( zs + idx_zs    );

		      rZZ0 = _mm256_mul_pd(rZZ0,rC);
		      rZZ0 = _mm256_mul_pd(rZZ0,rZS0);
		      rH0  = _mm256_add_pd(rH0,rZZ0);

		      _mm256_storeu_pd(grid + index_xy + k , rH0 );

		      idx_zs+=4; 
		      idx_zz+=4;
		    }
		}
	    }
	}
    }
}


// -----------------------------------------------------------------------------
static
void SE_grid_split_AVX_u8_gaussian_d(real* gmx_restrict grid, real* gmx_restrict q,
				     splinedata_t         * gmx_restrict spline, 
				     const pme_atomcomm_t * gmx_restrict atc,
				     const pmegrid_t      * gmx_restrict pmegrid)
{
  // unpack parameters
  const double* zs = (double*) spline->zs;
  const double* zx = (double*) spline->theta[0];
  const double* zy = (double*) spline->theta[1];
  const double* zz = (double*) spline->theta[2];
    
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

  __m256d rH0, rZZ0, rZS0, rC;
  __m256d rH1, rZZ1, rZS1;

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
    
    if(0){ // H[idx0] is 32-aligned
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	real qnzx = qn*zx[p*n+i];
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  rC = _mm256_set1_pd( qnzx*zy[p*n+j] );
	  idx_zz=p*n;
	  
	  for(k = 0; k<p; k+=8){
	    rH0  = _mm256_load_pd( grid + index_xy + k );
	    rH1  = _mm256_load_pd( grid + index_xy + k + 4 );
		      
	    rZZ0 = _mm256_load_pd( zz + idx_zz     );
	    rZZ1 = _mm256_load_pd( zz + idx_zz + 4 );
		      
	    rZS0 = _mm256_load_pd( zs + idx_zs    );
	    rZS1 = _mm256_load_pd( zs + idx_zs + 4);

	    rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
	    rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
		      
	    _mm256_store_pd( grid + index_xy + k    , rH0 );
	    _mm256_store_pd( grid + index_xy + k + 4, rH1 );
		      
	    idx_zs+=8; 
	    idx_zz+=8;
	  }
	}
      }
    }
    else{ // H[idx0] is 16-aligned, preventing nice vectorization
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	real qnzx = qn*zx[p*n+i];
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  rC = _mm256_set1_pd( qnzx*zy[p*n+j] );
	  idx_zz=p*n;
	  
	  for(k = 0; k<p; k+=8){
	    rH0  = _mm256_loadu_pd( grid + index_xy + k );
	    rH1  = _mm256_loadu_pd( grid + index_xy + k + 4 );
		      
	    rZZ0 = _mm256_load_pd( zz + idx_zz     );
	    rZZ1 = _mm256_load_pd( zz + idx_zz + 4 );
		      
	    rZS0 = _mm256_load_pd( zs + idx_zs    );
	    rZS1 = _mm256_load_pd( zs + idx_zs + 4);

	    rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
	    rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
		      
	    _mm256_storeu_pd( grid + index_xy + k    , rH0 );
	    _mm256_storeu_pd( grid + index_xy + k + 4, rH1 );
		      
	    idx_zs+=8; 
	    idx_zz+=8;
	  }
	}
      }
    }
  }
}

// -----------------------------------------------------------------------------
static
void SE_grid_split_AVX_P8_gaussian_d(real* gmx_restrict grid, real* gmx_restrict q,
				     splinedata_t         * gmx_restrict spline, 
				     const pme_atomcomm_t * gmx_restrict atc,
				     const pmegrid_t      * gmx_restrict pmegrid)
{
  // unpack parameters
  const double* zs = (double*) spline->zs;
  const double* zx = (double*) spline->theta[0];
  const double* zy = (double*) spline->theta[1];
  const double* zz = (double*) spline->theta[2];

  double qn;
  int idx_zs, i, j, n, nn;
  int index_x, index_xy;
  int i0,j0,k0;
  int * idxptr;

  int pny = pmegrid->s[YY];
  int pnz = pmegrid->s[ZZ];

  int offx = pmegrid->offset[XX];
  int offy = pmegrid->offset[YY];
  int offz = pmegrid->offset[ZZ];
   
  __m256d rZZ0, rZZ1; 
  __m256d rH0, rH1;
  __m256d rC, rZS0,rZS1;

  for(n=0; n<spline->n; n++){
    nn = spline->ind[n];
    qn = q[nn];
    idxptr = atc->idx[nn];
    i0 = idxptr[XX] - offx;
    j0 = idxptr[YY] - offy;
    k0 = idxptr[ZZ] - offz;
    idx_zs = 0;
    
    rZZ0 = _mm256_load_pd(zz + n*8     );
    rZZ1 = _mm256_load_pd(zz + n*8 + 4 );
    
    if(0){ // H[idx0] is 32-aligned
      for(i = 0; i<8; i++){
	index_x = (i0+i)*pny*pnz;
	real qnzx = qn*zx[8*n+i];
	for(j = 0; j<8; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  rC = _mm256_set1_pd( qnzx*zy[8*n+j] );

	  rH0  = _mm256_load_pd( grid + index_xy );
	  rH1  = _mm256_load_pd( grid + index_xy + 4);

	  rZS0 = _mm256_load_pd( zs + idx_zs);
	  rZS1 = _mm256_load_pd( zs + idx_zs + 4);

	  rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
	  rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
		  
	  _mm256_store_pd(grid + index_xy    ,  rH0);
	  _mm256_store_pd(grid + index_xy + 4,  rH1);

	  idx_zs += 8;
	}
      }
    }
    else{ // H[idx0] is 16-aligned, preventing nice vectorization
      for(i = 0; i<8; i++){
	index_x = (i0+i)*pny*pnz;
	real qnzx = qn*zx[8*n+i];
	for(j = 0; j<8; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  rC = _mm256_set1_pd( qnzx*zy[8*n+j] );

	  rH0  = _mm256_loadu_pd( grid + index_xy );
	  rH1  = _mm256_loadu_pd( grid + index_xy + 4);

	  rZS0 = _mm256_load_pd( zs + idx_zs);
	  rZS1 = _mm256_load_pd( zs + idx_zs + 4);

	  rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
	  rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
		  
	  _mm256_storeu_pd(grid + index_xy    ,  rH0);
	  _mm256_storeu_pd(grid + index_xy + 4,  rH1);

	  idx_zs += 8;
	}
      }      
    }
  }
}

// -------------------------------------------------------------------------
static
void SE_grid_split_AVX_P16_gaussian_d(real* gmx_restrict grid, real* gmx_restrict q,
				      splinedata_t         * gmx_restrict spline, 
				      const pme_atomcomm_t * gmx_restrict atc,
				      const pmegrid_t      * gmx_restrict pmegrid)
{
  // unpack parameters
  const double* zs = (double*) spline->zs;
  const double* zx = (double*) spline->theta[0];
  const double* zy = (double*) spline->theta[1];
  const double* zz = (double*) spline->theta[2];

  double qn;
  int idx0, idx_zs, i, j, n, nn;
  int index_x, index_xy;
  int i0,j0,k0;
  int * idxptr;

  int pny = pmegrid->s[YY];
  int pnz = pmegrid->s[ZZ];

  int offx = pmegrid->offset[XX];
  int offy = pmegrid->offset[YY];
  int offz = pmegrid->offset[ZZ];

  __m256d rZZ0, rZZ1, rZZ2, rZZ3; 
  __m256d rH0, rH1, rH2, rH3;
  __m256d rC, rZS0,rZS1,rZS2,rZS3;

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

    rZZ0 = _mm256_load_pd(zz + n*16     );
    rZZ1 = _mm256_load_pd(zz + n*16 + 4 );
    rZZ2 = _mm256_load_pd(zz + n*16 + 8 );
    rZZ3 = _mm256_load_pd(zz + n*16 + 12);
      
    if(0){ // H[idx0] is 32-aligned
      for(i = 0; i<16; i++){
	index_x = (i0+i)*pny*pnz;
	real qnzx = qn*zx[16*n+i];
	for(j = 0; j<16; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  rC = _mm256_set1_pd( qnzx*zy[16*n+j] );
		  
	  rH0  = _mm256_load_pd( grid + index_xy );
	  rH1  = _mm256_load_pd( grid + index_xy + 4);
	  rH2  = _mm256_load_pd( grid + index_xy + 8);
	  rH3  = _mm256_load_pd( grid + index_xy + 12);

	  rZS0 = _mm256_load_pd( zs + idx_zs);
	  rZS1 = _mm256_load_pd( zs + idx_zs + 4);
	  rZS2 = _mm256_load_pd( zs + idx_zs + 8);   
	  rZS3 = _mm256_load_pd( zs + idx_zs + 12);    

	  rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
	  rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
	  rH2 = _mm256_add_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rC),rZS2));
	  rH3 = _mm256_add_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rC),rZS3));

	  _mm256_store_pd(grid + index_xy    ,  rH0);
	  _mm256_store_pd(grid + index_xy + 4,  rH1);
	  _mm256_store_pd(grid + index_xy + 8,  rH2);
	  _mm256_store_pd(grid + index_xy + 12, rH3);

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
	  rC = _mm256_set1_pd( qnzx*zy[16*n+j] );
		  
	  rH0  = _mm256_loadu_pd( grid + index_xy );
	  rH1  = _mm256_loadu_pd( grid + index_xy + 4);
	  rH2  = _mm256_loadu_pd( grid + index_xy + 8);
	  rH3  = _mm256_loadu_pd( grid + index_xy + 12);

	  rZS0 = _mm256_load_pd( zs + idx_zs);
	  rZS1 = _mm256_load_pd( zs + idx_zs + 4);
	  rZS2 = _mm256_load_pd( zs + idx_zs + 8);   
	  rZS3 = _mm256_load_pd( zs + idx_zs + 12);    

	  rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
	  rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
	  rH2 = _mm256_add_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rC),rZS2));
	  rH3 = _mm256_add_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rC),rZS3));

	  _mm256_storeu_pd(grid + index_xy    ,  rH0);
	  _mm256_storeu_pd(grid + index_xy + 4,  rH1);
	  _mm256_storeu_pd(grid + index_xy + 8,  rH2);
	  _mm256_storeu_pd(grid + index_xy + 12, rH3);

	  idx_zs += 16;
	}
      }
    }
  }
}

/* =============================================================================
 * =============================== KAISER ROUTINES =============================
 * =============================================================================
 */
static
void SE_grid_split_AVX_P4_kaiser_d(real* gmx_restrict grid, real* gmx_restrict q,
				   splinedata_t         * gmx_restrict spline, 
				   const pme_atomcomm_t * gmx_restrict atc,
				   const pmegrid_t      * gmx_restrict pmegrid)
{
  // unpack parameters
  const double* zx = (double*) spline->theta[0];
  const double* zy = (double*) spline->theta[1];
  const double* zz = (double*) spline->theta[2];

  double qn;
  int i, j, nn;
  int index_x, index_xy;
  int i0,j0,k0;
  int * idxptr;

  int pny = pmegrid->s[YY];
  int pnz = pmegrid->s[ZZ];

  int offx = pmegrid->offset[XX];
  int offy = pmegrid->offset[YY];
  int offz = pmegrid->offset[ZZ];
  
  __m256d rZZ0;
  __m256d rH0;
  __m256d rC;

  for(int n=0; n<spline->n; n++)
    {
      nn = spline->ind[n];
      qn = q[nn];
      idxptr = atc->idx[nn];
      i0 = idxptr[XX] - offx;
      j0 = idxptr[YY] - offy;
      k0 = idxptr[ZZ] - offz;
    
      int n4 = n*4;
      rZZ0 = _mm256_load_pd(zz + n4);

      for(i = 0; i<4; i++)
	{
	  index_x = (i0+i)*pny*pnz;
	  double qnzx = qn*zx[n4+i];
	  for(j = 0; j<4; j++)
	    {
	      index_xy = index_x + (j0+j)*pnz + k0;
	      rC = _mm256_set1_pd( qnzx*zy[n4+j] );

	      rH0  = _mm256_loadu_pd( grid + index_xy );

#ifdef AVX_FMA
	      rH0 = _mm256_fmadd_pd(rZZ0,rC, rH0);
#else
	      rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(rZZ0,rC));
#endif
	      _mm256_storeu_pd( grid + index_xy , rH0);

	    }
	}
    }
}

static
void SE_grid_split_AVX_P6_kaiser_d(real* gmx_restrict grid, real* gmx_restrict q,
				   splinedata_t         * gmx_restrict spline, 
				   const pme_atomcomm_t * gmx_restrict atc,
				   const pmegrid_t      * gmx_restrict pmegrid)
{
  // unpack parameters
  const double* zx = (double*) spline->theta[0];
  const double* zy = (double*) spline->theta[1];
  const double* zz = (double*) spline->theta[2];

  double qn;
  int i, j, n;
  int index_x, index_xy;
  int i0,j0,k0;
  int * idxptr;

  int pny = pmegrid->s[YY];
  int pnz = pmegrid->s[ZZ];

  int offx = pmegrid->offset[XX];
  int offy = pmegrid->offset[YY];
  int offz = pmegrid->offset[ZZ];

  __m256d rZZ0; 
  __m256d rH0;
  __m256d rC0;

  __m128d rZZ1; 
  __m128d rH1;
  __m128d rC1;

  for(int nn=0; nn<spline->n; nn++)
    {
      n = spline->ind[nn];
      qn = q[n];
      idxptr = atc->idx[n];
      i0 = idxptr[XX] - offx;
      j0 = idxptr[YY] - offy;
      k0 = idxptr[ZZ] - offz;

      int n6 = n*6;
      
      rZZ0 = _mm256_loadu_pd(zz + n6     );
      rZZ1 = _mm_loadu_pd(  zz + n6 + 4 );

      for(i = 0; i<6; i++)
	{
	  index_x = (i0+i)*pny*pnz;
	  real qnzx = qn*zx[n6+i];
	  for(j = 0; j<6; j++)
	    {
	      index_xy = index_x + (j0+j)*pnz + k0;
	      rC0 = _mm256_set1_pd( qnzx*zy[n6+j] );
	      rC1 = _mm_set1_pd( qnzx*zy[n6+j] );

	      rH0  = _mm256_loadu_pd( grid + index_xy     );
	      rH1  = _mm_loadu_pd( grid + index_xy +  4);

#ifdef AVX_FMA
	      rH0 = _mm256_fmadd_pd(rZZ0,rC0, rH0);
	      rH1 = _mm_fmadd_pd(rZZ1,rC1, rH1);
#else
	      rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(rZZ0,rC0));
	      rH1 = _mm_add_pd(rH1,_mm_mul_pd(rZZ1,rC1));
#endif
	      _mm256_storeu_pd(grid + index_xy,      rH0);
	      _mm_storeu_pd(grid + index_xy + 4,  rH1);

	    }
	}
    }
}

// -----------------------------------------------------------------------------
static
void SE_grid_split_AVX_P8_kaiser_d(real* gmx_restrict grid, real* gmx_restrict q,
				   splinedata_t         * gmx_restrict spline, 
				   const pme_atomcomm_t * gmx_restrict atc,
				   const pmegrid_t      * gmx_restrict pmegrid)
{
  // unpack parameters
  const double* zx = (double*) spline->theta[0];
  const double* zy = (double*) spline->theta[1];
  const double* zz = (double*) spline->theta[2];

  double qn;
  int i, j, nn;
  int index_x, index_xy;
  int i0,j0,k0;
  int * idxptr;

  int pny = pmegrid->s[YY];
  int pnz = pmegrid->s[ZZ];

  int offx = pmegrid->offset[XX];
  int offy = pmegrid->offset[YY];
  int offz = pmegrid->offset[ZZ];

  __m256d rZZ0, rZZ1; 
  __m256d rH0, rH1;
  __m256d rC;

  for(int n=0; n<spline->n; n++)
    {
      nn = spline->ind[n];
      qn = q[nn];
      idxptr = atc->idx[nn];
      i0 = idxptr[XX] - offx;
      j0 = idxptr[YY] - offy;
      k0 = idxptr[ZZ] - offz;
      
      rZZ0 = _mm256_load_pd(zz + n*8     );
      rZZ1 = _mm256_load_pd(zz + n*8 + 4 );

      int n8 = n*8;
      for(i = 0; i<8; i++)
	{
	  index_x = (i0+i)*pny*pnz;
	  real qnzx = qn*zx[n8+i];
	  for(j = 0; j<8; j++)
	    {
	      index_xy = index_x + (j0+j)*pnz + k0;
	      rC = _mm256_set1_pd( qnzx*zy[n8+j] );

	      rH0  = _mm256_loadu_pd( grid + index_xy     );
	      rH1  = _mm256_loadu_pd( grid + index_xy +  4);

#ifdef AVX_FMA
	      rH0 = _mm256_fmadd_pd(rZZ0,rC, rH0);
	      rH1 = _mm256_fmadd_pd(rZZ1,rC, rH1);
#else
	      rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(rZZ0,rC));
	      rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(rZZ1,rC));
#endif
	      _mm256_storeu_pd(grid + index_xy,      rH0);
	      _mm256_storeu_pd(grid + index_xy + 4,  rH1);

	    }
	}
    }
}

// -----------------------------------------------------------------------------
static
void SE_grid_split_AVX_P16_kaiser_d(real* gmx_restrict grid, real* gmx_restrict q,
				    splinedata_t         * gmx_restrict spline, 
				    const pme_atomcomm_t * gmx_restrict atc,
				    const pmegrid_t      * gmx_restrict pmegrid)
{
  // unpack parameters
  const double* zx = (double*) spline->theta[0];
  const double* zy = (double*) spline->theta[1];
  const double* zz = (double*) spline->theta[2];

  double qn;
  int idx0, i, j, n;
  int index_x, index_xy;
  int i0,j0,k0;
  int * idxptr;

  int pny = pmegrid->s[YY];
  int pnz = pmegrid->s[ZZ];

  int offx = pmegrid->offset[XX];
  int offy = pmegrid->offset[YY];
  int offz = pmegrid->offset[ZZ];
  
  __m256d rZZ0, rZZ1, rZZ2, rZZ3; 
  __m256d rH0, rH1, rH2, rH3;
  __m256d rC;

  for(int nn=0; nn<spline->n; nn++)
    {
      n = spline->ind[nn];
      qn = q[n];
      idxptr = atc->idx[n];
      i0 = idxptr[XX] - offx;
      j0 = idxptr[YY] - offy;
      k0 = idxptr[ZZ] - offz;
      idx0 = i0*pny*pnz+j0*pnz+k0;
      _mm_prefetch( (void*) (grid+idx0), _MM_HINT_T0);
      
      int n16 = n*16;
	
      rZZ0 = _mm256_load_pd(zz + n16     );
      rZZ1 = _mm256_load_pd(zz + n16 + 4 );
      rZZ2 = _mm256_load_pd(zz + n16 + 8 );
      rZZ3 = _mm256_load_pd(zz + n16 + 12);

      for(i = 0 ; i<16; i++)
	{
	  index_x = (i0+i)*pny*pnz;
	  double qnzx = qn*zx[n16+i];
	  for(j = 0; j<16; j++)
	    {
	      index_xy = index_x + (j0+j)*pnz + k0;
	      rC = _mm256_set1_pd( qnzx*zy[n16+j] );

	      rH0  = _mm256_loadu_pd( grid + index_xy     );
	      rH1  = _mm256_loadu_pd( grid + index_xy +  4);
	      rH2  = _mm256_loadu_pd( grid + index_xy +  8);
	      rH3  = _mm256_loadu_pd( grid + index_xy + 12);

#ifdef AVX_FMA
	      rH0 = _mm256_fmadd_pd(rC,rZZ0,rH0);
	      rH1 = _mm256_fmadd_pd(rC,rZZ1,rH1);
	      rH2 = _mm256_fmadd_pd(rC,rZZ2,rH2);
	      rH3 = _mm256_fmadd_pd(rC,rZZ3,rH3);
#else
	      rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(rZZ0,rC));
	      rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(rZZ1,rC));
	      rH2 = _mm256_add_pd(rH2,_mm256_mul_pd(rZZ2,rC));
	      rH3 = _mm256_add_pd(rH3,_mm256_mul_pd(rZZ3,rC));
#endif
	      _mm256_storeu_pd(grid + index_xy,      rH0);
	      _mm256_storeu_pd(grid + index_xy + 4,  rH1);
	      _mm256_storeu_pd(grid + index_xy + 8,  rH2);
	      _mm256_storeu_pd(grid + index_xy + 12, rH3);

	    }
	}
    }
}

// ---------------------------------------------
static void
SE_grid_split_AVX_dispatch_d(real* grid, real* q, 
			     splinedata_t *spline,
			     const SE_params* params, 
			     const pme_atomcomm_t *atc,
			     const pmegrid_t *pmegrid,
			     const int se_set)
{
  const int p = params->P;
  const int incrj = params->dims[2]; // middle increment
  const int incri = params->npdims[2]*(params->dims[1]);// outer increment
    
  if(se_set==1)
    {
      // if P, or either increments are not divisible by 4,fall back on vanilla
      if( isnot_div_by_4(p) || isnot_div_by_4(incri) || isnot_div_by_4(incrj) ) {
	__DISPATCHER_MSG("[FGG GRID AVX DOUBLE] AVX Abort (PARAMS)\n");
	SE_grid_split_SSE_dispatch_d(grid, q, spline, params, atc, pmegrid, se_set);
	return;
      }
    }
  else
    {
      // if either increments are not divisible by 4,fall back on vanilla
      if( isnot_div_by_4(incri) || isnot_div_by_4(incrj) ) {
	__DISPATCHER_MSG("[FKG GRID AVX DOUBLE] AVX Abort (PARAMS)\n");
      SE_grid_split_SSE_dispatch_d(grid, q, spline, params, atc, pmegrid, se_set);
      return;
      }
    }
  
  if(se_set==1)
    {
      // otherwise the preconditions for AVX codes are satisfied. 
      if(p==16){
	// specific for p=16
	__DISPATCHER_MSG("[FGG GRID AVX DOUBLE] P=16\n");
	SE_grid_split_AVX_P16_gaussian_d(grid, q, spline, atc, pmegrid); 
      }
      else if(p==8){
	// specific for p=8
	__DISPATCHER_MSG("[FGG GRID AVX DOUBLE] P=8\n");
	SE_grid_split_AVX_P8_gaussian_d(grid, q, spline, atc, pmegrid); 
      }
      else if(p%8==0){
	// specific for p divisible by 8
	__DISPATCHER_MSG("[FGG GRID AVX DOUBLE] P unroll 8\n");
	SE_grid_split_AVX_u8_gaussian_d(grid, q, spline, atc, pmegrid); 
      }
      else if(p%4==0){
	// specific for p divisible by 4
	__DISPATCHER_MSG("[FGG GRID AVX DOUBLE] P unroll 4\n");
	SE_grid_split_AVX_gaussian_d(grid, q, spline, atc, pmegrid); 
      }
      else{
	// vanilla SSE code (any even p)
	__DISPATCHER_MSG("[FGG GRID AVX DOUBLE] Vanilla\n");
	SE_grid_split_SSE_gaussian_d(grid, q, spline, atc, pmegrid);
      }
    } // se_set==1 Gaussian
  else // Kaiser
    {
      // otherwise the preconditions for AVX codes are satisfied.
      if(p==16){
      	// specific for p=16
      	__DISPATCHER_MSG("[FKG GRID AVX DOUBLE] P=16\n");
      	SE_grid_split_AVX_P16_kaiser_d(grid, q, spline, atc, pmegrid);
      }
      else if(p==8){
      	// specific for p=8
      	__DISPATCHER_MSG("[FKG GRID AVX DOUBLE] P=8\n");
      	SE_grid_split_AVX_P8_kaiser_d(grid, q, spline, atc, pmegrid);
      }
      else if(p==6){
      	// specific for p=6
      	__DISPATCHER_MSG("[FKG GRID AVX DOUBLE] P=6\n");
      	SE_grid_split_AVX_P6_kaiser_d(grid, q, spline, atc, pmegrid);
      }      
      else if(p==4){
      	// specific for p=4
      	__DISPATCHER_MSG("[FKG GRID AVX DOUBLE] P unroll 4\n");
      	SE_grid_split_AVX_P4_kaiser_d(grid, q, spline, atc, pmegrid);
      }
      else{
      	// vanilla SSE code (any even p)
      	__DISPATCHER_MSG("[FKG GRID AVX DOUBLE] AVX Abort\n");
      	SE_grid_split_SSE_kaiser_d(grid, q, spline, atc, pmegrid);
      }
    }
}


#endif //GMX_SIMD_X86_AVX_256
#endif //GMX_DOUBLE
#endif //_SE_GRID_AVX_256_DOUBLE_
