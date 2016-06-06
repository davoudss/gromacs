#ifndef _SE_INT_SSE_SINGLE_
#define _SE_INT_SSE_SINGLE_

/* SE SSE single integration */
#if GMX_DOUBLE==0
#include "se.h"

// -----------------------------------------------------------------------------
static void 
SE_int_split(rvec* force,  real* grid, real* q,
	     splinedata_t *spline,
	     const SE_FGG_params* params, real scale, 
	     gmx_bool bClearF,
	     const pme_atomcomm_t *atc,
	     const gmx_pme_t *pme)
{
  // unpack params
  const float*   H   = (float*) grid;
  const float*   zs  = (float*) spline->zs;
  const float*   zx  = (float*) spline->theta[0];
  const float*   zy  = (float*) spline->theta[1];
  const float*   zz  = (float*) spline->theta[2];
  const float*   zfx = (float*) spline->dtheta[0];
  const float*   zfy = (float*) spline->dtheta[1];
  const float*   zfz = (float*) spline->dtheta[2];


  const int   p = params->P;
  const float h = params->h;

  int i,j,k,m,idx_zs,idx_zz,mm;
  int * idxptr;
  int i0,j0,k0;
  int index_x, index_xy;
  float force_m[3], cij, Hzc, qm, h3 = h*h*h;
#ifdef CALC_ENERGY
  float phi_m;
#endif

  int pny   = pme->pmegrid_ny;
  int pnz   = pme->pmegrid_nz;

  for(m=0; m<spline->n; m++)
    {
      mm  = spline->ind[m];
      qm = q[mm];
      idxptr = atc->idx[mm];
      i0 = idxptr[XX];
      j0 = idxptr[YY];
      k0 = idxptr[ZZ];

      force_m[0] = 0; force_m[1] = 0; force_m[2] = 0;
#ifdef CALC_ENERGY
      phi_m = 0;
#endif
      idx_zs = 0;

      for(i = 0; i<p; i++)
	{
	  index_x = (i0+i)*pny*pnz;
	  for(j = 0; j<p; j++)
	    {
	      cij = zx[m*p+i]*zy[m*p+j];
	      idx_zz=m*p;
	      index_xy = index_x + (j0+j)*pnz;
	      for(k = 0; k<p; k++)
		{
		  Hzc = H[index_xy+(k0+k)]*zs[idx_zs]*zz[idx_zz]*cij;
#ifdef CALC_ENERGY
		  phi_m      += Hzc;
#endif
		  force_m[0] += Hzc*zfx[m*p+i];
		  force_m[1] += Hzc*zfy[m*p+j];
		  force_m[2] += Hzc*zfz[m*p+k];
		    
		  idx_zs++; idx_zz++;
		}
	    }
	}

      if(bClearF)
      	{
      	  force[m][XX] = 0;
      	  force[m][YY] = 0;
      	  force[m][ZZ] = 0;
      	}

      force[m][XX] += -qm*scale*h3*force_m[0];
      force[m][YY] += -qm*scale*h3*force_m[1];
      force[m][ZZ] += -qm*scale*h3*force_m[2];

#ifdef CALC_ENERGY
      st->phi[m]   = -h3*scale*phi_m;
#endif
    }
}

// -----------------------------------------------------------------------------
static void 
SE_int_split_SSE(rvec *force, real *grid, real *q,
		 splinedata_t *spline,
		 const SE_FGG_params *params, real scale,
		 gmx_bool bClearF,
		 const pme_atomcomm_t *atc,
		 const gmx_pme_t *pme)
{

  // unpack params
  const float*    H = (float*) grid;
  const float*   zs = (float*) spline->zs;
  const float*   zx = (float*) spline->theta[0];
  const float*   zy = (float*) spline->theta[1];
  const float*   zz = (float*) spline->theta[2];
  const float*  zfx = (float*) spline->dtheta[0];
  const float*  zfy = (float*) spline->dtheta[1];
  const float*  zfz = (float*) spline->dtheta[2];

  const int   p = params->P;
  const int   N = params->N;
  const float h = params->h;

  int i,j,k,m,idx,idx_zs,idx_zz, mm;
  int * idxptr;
  int i0,j0,k0;
  int index_x, index_xy;

  float qm, h3=h*h*h;
  float sx[4] MEM_ALIGNED;
  float sy[4] MEM_ALIGNED;
  float sz[4] MEM_ALIGNED;

  int pny   = pme->pmegrid_ny;
  int pnz   = pme->pmegrid_nz;

  __m128 rH0, rZZ0, rZS0,rZFZ0;
  __m128 rC, rCX, rCY;
  __m128 rFX, rFY, rFZ;

#ifdef CALC_ENERGY
  float s[4] MEM_ALIGNED;
  __m128 rP, rCP;
#endif

  for(m=0; m<N; m++){
    mm = spline->ind[m];
    qm = q[mm];
    idxptr = atc->idx[mm];
    i0 = idxptr[XX];
    j0 = idxptr[YY];
    k0 = idxptr[ZZ];
    idx = i0*pny*pnz+j0*pnz+k0;
    
    idx_zs = 0;
    rFX = _mm_setzero_ps();
    rFY = _mm_setzero_ps();
    rFZ = _mm_setzero_ps();
#ifdef CALC_ENERGY
    rP  = _mm_setzero_ps();
#endif

    if(idx%4==0){ // H[idx] is 16-aligned so vectorization simple
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz;
	  rC  = _mm_set1_ps( zx[m*p+i]*zy[m*p+j]);
	  rCX = _mm_set1_ps( zx[m*p+i]*zy[m*p+j]*zfx[m*p+i]);
	  rCY = _mm_set1_ps( zx[m*p+i]*zy[m*p+j]*zfy[m*p+j]);
#ifdef CALC_ENERGY
	  rCP = _mm_set1_ps( zx[m*p+i]*zy[m*p+j]);
#endif
	  
	  idx_zz=m*p;
	  for(k = 0; k<p; k+=4){
	    rZFZ0= _mm_load_ps( zfz+ m*p+k );
	    rH0  = _mm_load_ps( H  + index_xy + k0 + k );
	    rZZ0 = _mm_load_ps( zz + idx_zz);
	    rZS0 = _mm_load_ps( zs + idx_zs);
	    rFX = _mm_add_ps(rFX,_mm_mul_ps(rH0,_mm_mul_ps(_mm_mul_ps(rZZ0,rCX),rZS0)));
	    rFY = _mm_add_ps(rFY,_mm_mul_ps(rH0,_mm_mul_ps(_mm_mul_ps(rZZ0,rCY),rZS0)));
	    rFZ = _mm_add_ps(rFZ,_mm_mul_ps(_mm_mul_ps(rH0,_mm_mul_ps(_mm_mul_ps(rZZ0,rC),rZS0)),rZFZ0));

#ifdef CALC_ENERGY
	    rP = _mm_add_ps(rP,_mm_mul_ps(_mm_mul_ps(rH0,_mm_mul_ps(_mm_mul_ps(rZZ0,rCP),rZS0)),rZFZ0));
#endif
			
	    idx_zs+=4; 
	    idx_zz+=4;
	  }
	}
      }
    }
    else{ // H[idx] not 16-aligned, so use non-aligned loads
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz;
	  rC  = _mm_set1_ps( zx[m*p+i]*zy[m*p+j] );
	  rCX = _mm_set1_ps( zx[m*p+i]*zy[m*p+j]*zfx[m*p+i] );
	  rCY = _mm_set1_ps( zx[m*p+i]*zy[m*p+j]*zfy[m*p+j] );
#ifdef CALC_ENERGY
	  rCP = _mm_set1_ps( zx[m*p+i]*zy[m*p+j]);
#endif

	  idx_zz=m*p;
	  for(k = 0; k<p; k+=4){
	    rZFZ0= _mm_load_ps( zfz + m*p+k );
	    rH0  = _mm_loadu_ps( H+index_xy + k0 + k );
	    rZZ0 = _mm_load_ps( zz + idx_zz);
	    rZS0 = _mm_load_ps( zs + idx_zs);
	    rFX = _mm_add_ps(rFX,_mm_mul_ps(rH0,_mm_mul_ps(_mm_mul_ps(rZZ0,rCX),rZS0)));
	    rFY = _mm_add_ps(rFY,_mm_mul_ps(rH0,_mm_mul_ps(_mm_mul_ps(rZZ0,rCY),rZS0)));
	    rFZ = _mm_add_ps(rFZ,_mm_mul_ps(_mm_mul_ps(rH0,_mm_mul_ps(_mm_mul_ps(rZZ0,rC),rZS0)),rZFZ0));

#ifdef CALC_ENERGY
	    rP = _mm_add_ps(rP,_mm_mul_ps(_mm_mul_ps(rH0,_mm_mul_ps(_mm_mul_ps(rZZ0,rCP),rZS0)),rZFZ0));
#endif
			
	    idx_zs+=4; 
	    idx_zz+=4;
	  }
	}
      }

    }
    _mm_store_ps(sx,rFX);
    _mm_store_ps(sy,rFY);
    _mm_store_ps(sz,rFZ);

    if(bClearF){
      force[m][XX] = 0;
      force[m][YY] = 0;
      force[m][ZZ] = 0;
    }

    force[m][XX] += -qm*scale*h3*(sx[0]+sx[1]+sx[2]+sx[3]);
    force[m][YY] += -qm*scale*h3*(sy[0]+sy[1]+sy[2]+sy[3]);
    force[m][ZZ] += -qm*scale*h3*(sz[0]+sz[1]+sz[2]+sz[3]);

#ifdef CALC_ENERGY
    _mm_store_ps(s,rP);
    st->phi[m] = -scale*h3*(s[0]+s[1]+s[2]+s[3]);
#endif
  }
}

// -----------------------------------------------------------------------------
static void
SE_int_split_SSE_P8(rvec *force, real *grid, real *q,
		    splinedata_t *spline,
		    const SE_FGG_params *params, real scale,
		    gmx_bool bClearF,
		    const pme_atomcomm_t *atc,
		    const gmx_pme_t *pme)
{
  // unpack params
  const float*    H = (float*) grid;
  const float*   zs = (float*) spline->zs;
  const float*   zx = (float*) spline->theta[0];
  const float*   zy = (float*) spline->theta[1];
  const float*   zz = (float*) spline->theta[2];
  const float*  zfx = (float*) spline->dtheta[0];
  const float*  zfy = (float*) spline->dtheta[1];
  const float*  zfz = (float*) spline->dtheta[2];

  /* ASSUME P=8 const int p = params->P; */
  const int   N = params->N;
  const float h = params->h;

  int i,j,m,idx,idx_zs, mm;
  int * idxptr;
  int i0,j0,k0;
  int index_x, index_xy;

  float qm, h3=h*h*h;
  float sx[4] MEM_ALIGNED;
  float sy[4] MEM_ALIGNED;
  float sz[4] MEM_ALIGNED;

  int pny   = pme->pmegrid_ny;
  int pnz   = pme->pmegrid_nz;

  // hold entire zz vector
  __m128 rZZ0, rZZ1; 
  __m128 rC, rCX, rCY;
  __m128 rH0, rH1; 
  __m128 rZS0, rZS1;
  __m128 rFX, rFY, rFZ,rZFZ0,rZFZ1;
#ifdef CALC_ENERGY
  float s[4] MEM_ALIGNED;
  __m128 rP, rCP;
#endif

  for(m=0; m<N; m++){
    mm = spline->ind[m];
    qm = q[mm];

    idxptr = atc->idx[mm];
    i0 = idxptr[XX];
    j0 = idxptr[YY];
    k0 = idxptr[ZZ];
    idx = i0*pny*pnz+j0*pnz+k0;
    
    idx_zs = 0;
    rFX = _mm_setzero_ps();
    rFY = _mm_setzero_ps();
    rFZ = _mm_setzero_ps();
#ifdef CALC_ENERGY
    rP  = _mm_setzero_ps();
#endif
    
    /* hoist load of ZZ vector */
    rZZ0  = _mm_load_ps(zz  + m*8     );
    rZZ1  = _mm_load_ps(zz  + m*8 + 4 );
    rZFZ0 = _mm_load_ps(zfz + m*8     );
    rZFZ1 = _mm_load_ps(zfz + m*8 + 4 );

    if(idx%4==0){ // H[idx] is 16-aligned so vectorization simple
      for(i = 0; i<8; i++){
	index_x = (i0+i)*pny*pnz;
	for(j = 0; j<8; j++){
	  index_xy = index_x + (j0+j)*pnz;
	  rC = _mm_set1_ps( zx[m*8+i]*zy[m*8+j]);
	  rCX = _mm_set1_ps(zx[m*8+i]*zy[m*8+j]*zfx[m*8+i]);
	  rCY = _mm_set1_ps(zx[m*8+i]*zy[m*8+j]*zfy[m*8+j]);
#ifdef CALC_ENERGY
	  rCP= _mm_set1_ps( zx[m*8+i]*zy[m*8+j]);
#endif

	  rH0  = _mm_load_ps( H+index_xy + k0    );
	  rH1  = _mm_load_ps( H+index_xy + k0 + 4);

	  rZS0 = _mm_load_ps( zs + idx_zs    );
	  rZS1 = _mm_load_ps( zs + idx_zs + 4);
		   
	  rFX = _mm_add_ps(rFX,_mm_mul_ps(rH0,_mm_mul_ps(_mm_mul_ps(rZZ0,rCX),rZS0)));
	  rFX = _mm_add_ps(rFX,_mm_mul_ps(rH1,_mm_mul_ps(_mm_mul_ps(rZZ1,rCX),rZS1)));
	  rFY = _mm_add_ps(rFY,_mm_mul_ps(rH0,_mm_mul_ps(_mm_mul_ps(rZZ0,rCY),rZS0)));
	  rFY = _mm_add_ps(rFY,_mm_mul_ps(rH1,_mm_mul_ps(_mm_mul_ps(rZZ1,rCY),rZS1)));
	  rFZ = _mm_add_ps(rFZ,_mm_mul_ps(rH0,_mm_mul_ps(_mm_mul_ps(_mm_mul_ps(rZFZ0,rZZ0),rC),rZS0)));
	  rFZ = _mm_add_ps(rFZ,_mm_mul_ps(rH1,_mm_mul_ps(_mm_mul_ps(_mm_mul_ps(rZFZ1,rZZ1),rC),rZS1)));
#ifdef CALC_ENERGY
	  rP  = _mm_add_ps(rP,_mm_mul_ps(rH0,_mm_mul_ps(_mm_mul_ps(_mm_mul_ps(rZFZ0,rZZ0),rCP),rZS0)));
	  rP = _mm_add_ps(rP,_mm_mul_ps(rH1,_mm_mul_ps(_mm_mul_ps(_mm_mul_ps(rZFZ1,rZZ1),rCP),rZS1)));
#endif

	  idx_zs +=8;
	}
      }
    }
    else{ // H[idx] not 16-aligned, so use non-aligned loads
      for(i = 0; i<8; i++){
	index_x = (i0+i)*pny*pnz;
	for(j = 0; j<8; j++){
	  index_xy = index_x + (j0+j)*pnz;
	  rC = _mm_set1_ps( zx[m*8+i]*zy[m*8+j]);
	  rCX = _mm_set1_ps(zx[m*8+i]*zy[m*8+j]*zfx[m*8+i]);
	  rCY = _mm_set1_ps(zx[m*8+i]*zy[m*8+j]*zfy[m*8+j]);
#ifdef CALC_ENERGY
	  rCP= _mm_set1_ps( zx[m*8+i]*zy[m*8+j]);
#endif

	  rH0  = _mm_loadu_ps( H+index_xy + k0    );
	  rH1  = _mm_loadu_ps( H+index_xy + k0 + 4);
		
	  rZS0 = _mm_load_ps( zs + idx_zs    );
	  rZS1 = _mm_load_ps( zs + idx_zs + 4);
		    
	  rFX = _mm_add_ps(rFX,_mm_mul_ps(rH0,_mm_mul_ps(_mm_mul_ps(rZZ0,rCX),rZS0)));
	  rFX = _mm_add_ps(rFX,_mm_mul_ps(rH1,_mm_mul_ps(_mm_mul_ps(rZZ1,rCX),rZS1)));
	  rFY = _mm_add_ps(rFY,_mm_mul_ps(rH0,_mm_mul_ps(_mm_mul_ps(rZZ0,rCY),rZS0)));
	  rFY = _mm_add_ps(rFY,_mm_mul_ps(rH1,_mm_mul_ps(_mm_mul_ps(rZZ1,rCY),rZS1)));
	  rFZ = _mm_add_ps(rFZ,_mm_mul_ps(rH0,_mm_mul_ps(_mm_mul_ps(_mm_mul_ps(rZFZ0,rZZ0),rC),rZS0)));
	  rFZ = _mm_add_ps(rFZ,_mm_mul_ps(rH1,_mm_mul_ps(_mm_mul_ps(_mm_mul_ps(rZFZ1,rZZ1),rC),rZS1)));
#ifdef CALC_ENERGY
	  rP  = _mm_add_ps(rP,_mm_mul_ps(rH0,_mm_mul_ps(_mm_mul_ps(_mm_mul_ps(rZFZ0,rZZ0),rCP),rZS0)));
	  rP = _mm_add_ps(rP,_mm_mul_ps(rH1,_mm_mul_ps(_mm_mul_ps(_mm_mul_ps(rZFZ1,rZZ1),rCP),rZS1)));
#endif

	  idx_zs +=8;
	}
      }
    }

    _mm_store_ps(sx,rFX);
    _mm_store_ps(sy,rFY);
    _mm_store_ps(sz,rFZ);

    if(bClearF){
      force[m][XX] = 0;
      force[m][YY] = 0;
      force[m][ZZ] = 0;
    }

    force[m][XX] += -qm*scale*h3*(sx[0]+sx[1]+sx[2]+sx[3]);
    force[m][YY] += -qm*scale*h3*(sy[0]+sy[1]+sy[2]+sy[3]);
    force[m][ZZ] += -qm*scale*h3*(sz[0]+sz[1]+sz[2]+sz[3]);

#ifdef CALC_ENERGY
    _mm_store_ps(s,rP);
    st->phi[m] = -scale*h3*(s[0]+s[1]+s[2]+s[3]);
#endif
  }
}

// -----------------------------------------------------------------------------
static void 
SE_int_split_SSE_dispatch(rvec *force, real *grid, real *q,
			  splinedata_t *spline,
			  const SE_FGG_params *params, real scale,
			  gmx_bool bClearF,
			  const pme_atomcomm_t *atc,
			  const gmx_pme_t *pme
			  )
{
  const int p = params->P;
  const int incrj = params->dims[2]; // middle increment
  const int incri = params->npdims[2]*(params->dims[1]);// outer increment

  // if P is odd, or if either increment is odd, fall back on vanilla
  if( isnot_div_by_4(p) || isnot_div_by_4(incri) || isnot_div_by_4(incrj)  || (p%4)!=0 )
    {
      __DISPATCHER_MSG("[FGG INT SSE SINGLE] SSE Abort (PARAMS)\n");
      SE_int_split(force, grid, q, spline, params, scale, bClearF,
		   atc,pme);
      return;
    }
    
  // otherwise the preconditions for SSE codes are satisfied. 
  if(p==8){
    // specific for p=8
    __DISPATCHER_MSG("[FGG INT SSE SINGLE] P=8\n");
    SE_int_split_SSE_P8(force, grid, q, spline, params, scale, bClearF,
			atc, pme);
  } 
  else{
    // vanilla SSE code for p divisible by 4
    __DISPATCHER_MSG("[FGG INT SSE SINGLE] Vanilla\n");
    SE_int_split_SSE(force, grid, q, spline, params, scale, bClearF,
		     atc, pme);
  }
}

#endif // not GMX_DOUBLE
#endif //_SE_INT_SSE_SINGLE_
