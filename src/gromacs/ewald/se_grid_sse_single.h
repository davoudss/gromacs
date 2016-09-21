#ifndef _SE_GRID_SSE_SINGLE_
#define _SE_GRID_SSE_SINGLE_

/* SE int SSE single gridding */
#if GMX_DOUBLE == 0
#include "se.h"

// -----------------------------------------------------------------------------
static void SE_grid_split(real* grid, real* q,
			  splinedata_t *spline,
			  const pme_atomcomm_t *atc,
			  const pmegrid_t* pmegrid)
{
  // unpack parameters
  const float*   zs = (float*) spline->zs;
  const float*   zx = (float*) spline->theta[0];
  const float*   zy = (float*) spline->theta[1];
  const float*   zz = (float*) spline->theta[2];
    
  const int p = pmegrid->order-1;

  float cij0,qn;
  int zidx, idxzz, i, j, k, n, nn;
  int index_x, index_xy,index_xyz;
  int i0,j0,k0;
  int * idxptr;

  int pny = pmegrid->s[YY];
  int pnz = pmegrid->s[ZZ];

  int offx = pmegrid->offset[XX];
  int offy = pmegrid->offset[YY];
  int offz = pmegrid->offset[ZZ];

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
	  for(j = 0; j<p; j++)
	    {
	      cij0 = zx[p*n+i]*zy[p*n+j];
	      idxzz=p*n;
	      index_xy = index_x + (j0+j)*pnz;
	      for(k = 0; k<p; k++)
		{
		  index_xyz = index_xy + (k0+k);
		  grid[index_xyz] += zs[zidx]*zz[idxzz]*cij0*qn;
		  zidx++; idxzz++;
		}
	    }
	}
    }
}


// -----------------------------------------------------------------------------
static void SE_grid_split_SSE(real* grid, real* q,
			      splinedata_t *spline,
			      const pme_atomcomm_t *atc,
			      const pmegrid_t* pmegrid)
{
  // unpack parameters
  const float*   zs = (float*) spline->zs;
  const float*   zx = (float*) spline->theta[0];
  const float*   zy = (float*) spline->theta[1];
  const float*   zz = (float*) spline->theta[2];
    
  const int p = pmegrid->order-1;

  int idx0, idx_zs, idx_zz, i, j, k, n, nn;
  int index_x, index_xy;
  int i0,j0,k0;
  int * idxptr;
  float qn;

  int pny = pmegrid->s[YY];
  int pnz = pmegrid->s[ZZ];

  int offx = pmegrid->offset[XX];
  int offy = pmegrid->offset[YY];
  int offz = pmegrid->offset[ZZ];

  __m128 rH0, rZZ0, rZS0, rC;

  for(n=0; n<spline->n; n++){
    nn = spline->ind[n];
    qn = q[nn];
    idxptr = atc->idx[nn];
    i0 = idxptr[XX] - offx;
    j0 = idxptr[YY] - offy;
    k0 = idxptr[ZZ] - offz;
    idx0 = i0*pny*pnz+j0*pnz+k0;

    idx_zs = 0;
    
    if(idx0%4 == 0){ // H[idx0] is 16-aligned
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz;
	  rC = _mm_set1_ps( qn*zx[p*n+i]*zy[p*n+j] );
	  idx_zz=p*n;
	  for(k = 0; k<p; k+=4){
	    rH0  = _mm_load_ps( grid+index_xy + k0 + k  );
	    rZZ0 = _mm_load_ps( zz + idx_zz     );
	    rZS0 = _mm_load_ps( zs + idx_zs    );
	    
	    rZZ0 = _mm_mul_ps(rZZ0,rC);
	    rZZ0 = _mm_mul_ps(rZZ0,rZS0);
	    rH0  = _mm_add_ps(rH0,rZZ0);
	    
	    _mm_store_ps( grid +index_xy + k0 + k , rH0 );
	    
	    idx_zs+=4; 
	    idx_zz+=4;
	  }
	}
      }
    }
    else{ // H[idx0] is 8-aligned, preventing nice vectorization
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz;
	  rC = _mm_set1_ps( qn*zx[p*n+i]*zy[p*n+j] );
	  idx_zz=p*n;
	  for(k = 0; k<p; k+=4){
	    rH0  = _mm_loadu_ps( grid+index_xy + k0 + k );
	    rZZ0 = _mm_load_ps( zz + idx_zz );
	    rZS0 = _mm_load_ps( zs + idx_zs );
	    rZZ0 = _mm_mul_ps(rZZ0,rC);
	    rZZ0 = _mm_mul_ps(rZZ0,rZS0);
	    rH0  = _mm_add_ps(rH0,rZZ0);
	    _mm_storeu_ps( grid+index_xy + k0 + k, rH0 );
	    
	    idx_zs+=4;
	    idx_zz+=4;
	  }
	}
      }
    }
  }
}

// -----------------------------------------------------------------------------
static void SE_grid_split_SSE_P8(real *grid, real *q,
				 splinedata_t *spline,
				 const pme_atomcomm_t *atc,
				 const pmegrid_t* pmegrid)
{
  // unpack parameters
  const float*   zs = (float*) spline->zs;
  const float*   zx = (float*) spline->theta[0];
  const float*   zy = (float*) spline->theta[1];
  const float*   zz = (float*) spline->theta[2];

  float qn;
  int idx, idx_zs, i, j, n , nn;

  int index_x, index_xy;
  int i0,j0,k0;
  int * idxptr;
  int pny = pmegrid->s[YY];
  int pnz = pmegrid->s[ZZ];

  int offx = pmegrid->offset[XX];
  int offy = pmegrid->offset[YY];
  int offz = pmegrid->offset[ZZ];

  __m128 rZZ0, rZZ1; 
  __m128 rH0, rH1;
  __m128 rC, rZS0;

  for(n=0; n<spline->n; n++){
    nn = spline->ind[n];
    qn = q[nn];
    idxptr = atc->idx[nn];
    i0 = idxptr[XX] - offx;
    j0 = idxptr[YY] - offy;
    k0 = idxptr[ZZ] - offz;
    idx = i0*pny*pnz+j0*pnz+k0;

    _mm_prefetch( (void*) (H+idx), _MM_HINT_T0);
    
    idx_zs = 0;
    _mm_prefetch( (void*) zs, _MM_HINT_T0);
    
    rZZ0 = _mm_load_ps(zz + n*8     );
    rZZ1 = _mm_load_ps(zz + n*8 + 4 );
    if(idx%4 == 0){ // H[idx0] is 16-aligned
      for(i = 0; i<8; i++){
	index_x = (i0+i)*pny*pnz;
	for(j = 0; j<8; j++){
	  index_xy = index_x + (j0+j)*pnz;

	  rC = _mm_set1_ps( qn*zx[8*n+i]*zy[8*n+j] );
	  
	  rH0  = _mm_load_ps( grid+index_xy + k0    );
	  rH1  = _mm_load_ps( grid+index_xy + k0 + 4);
	  rZS0 = _mm_load_ps( zs + idx_zs);
	  rH0 = _mm_add_ps(rH0,_mm_mul_ps(_mm_mul_ps(rZZ0,rC),rZS0));
	  
	  rZS0 = _mm_load_ps( zs + idx_zs + 4);                   
	  rH1 = _mm_add_ps(rH1,_mm_mul_ps(_mm_mul_ps(rZZ1,rC),rZS0));
	  
	  _mm_store_ps(grid + index_xy + k0    , rH0);
	  _mm_store_ps(grid + index_xy + k0 + 4, rH1);
	  
	  idx_zs += 8;
	}
      }
    }
    else{ // H[idx0] is 8-aligned, preventing nice vectorization
      for(i = 0; i<8; i++){
	index_x = (i0+i)*pny*pnz;
	for(j = 0; j<8; j++){
	  index_xy = index_x + (j0+j)*pnz;
	  rC = _mm_set1_ps( qn*zx[8*n+i]*zy[8*n+j] );
	  
	  rH0  = _mm_loadu_ps( grid+index_xy + k0    );
	  rH1  = _mm_loadu_ps( grid+index_xy + k0 + 4);
	  
	  // if zs does not have 16-byte alignment, this will core.
	  // PLATFORM AND COMPILER DEPENDENT (FIXME)
	  rZS0 = _mm_load_ps( zs + idx_zs);
	  rH0 = _mm_add_ps(rH0,_mm_mul_ps(_mm_mul_ps(rZZ0,rC),rZS0));
	  
	  rZS0 = _mm_load_ps( zs + idx_zs + 4);                   
	  rH1 = _mm_add_ps(rH1,_mm_mul_ps(_mm_mul_ps(rZZ1,rC),rZS0));
	  
	  _mm_storeu_ps(grid + index_xy + k0     , rH0);
	  _mm_storeu_ps(grid + index_xy + k0  + 4, rH1);
	  
	  idx_zs += 8;
	}
      }
    }
  }
}

// -----------------------------------------------------------------------------
static void 
SE_grid_split_SSE_dispatch(real* grid, real* q, 
			   splinedata_t *spline,
			   const SE_FGG_params* params,
			   const pme_atomcomm_t *atc,
			   const pmegrid_t *pmegrid)
{
  const int p = params->P;
  const int incrj = params->dims[2]; // middle increment
  const int incri = params->npdims[2]*(params->dims[1]);// outer increment

  // if P is odd, or if either increment is odd, fall back on vanilla
  if( isnot_div_by_4(p) || isnot_div_by_4(incri) || isnot_div_by_4(incrj)  || (p%4)!=0)
    {
      __DISPATCHER_MSG("[FGG GRID SSE SINGLE] SSE Abort (PARAMS)\n");
      SE_grid_split(grid, q, spline, atc,pmegrid);
      return;
    }
   
  // otherwise the preconditions for SSE codes are satisfied.    
  if(p==8){
    // specific for p=8
    __DISPATCHER_MSG("[FGG GRID SSE SINGLE] P=8\n");
    SE_grid_split_SSE_P8(grid, q, spline, atc, pmegrid);
  }
  else {
    // vanilla SSE code for p divisible by 4
    __DISPATCHER_MSG("[FGG GRID SSE SINGLE] Vanilla\n");
    SE_grid_split_SSE(grid, q, spline, atc, pmegrid);
  }
}

#endif // not GMX_DOUBLE

#endif //_SE_GRID_SSE_SINGLE_
