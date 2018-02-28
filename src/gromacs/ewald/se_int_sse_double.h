#ifndef _SE_INT_SSE_DOUBLE_
#define _SE_INT_SSE_DOUBLE_

/* SE int SSE double integration */
#if GMX_DOUBLE==1
#include "se.h"

static inline double
reduce_sse(__m128d a)
{
  a = _mm_hadd_pd(a,a);
  return *reinterpret_cast<double *>(&a);
}

static void 
SE_int_split_d(rvec * gmx_restrict force,  real * gmx_restrict grid, 
	       real * gmx_restrict q,
	       splinedata_t         * gmx_restrict spline,
	       const SE_params  * gmx_restrict params, 
	       real scale, gmx_bool bClearF,
	       const pme_atomcomm_t * gmx_restrict atc,
	       const gmx_pme_t      * gmx_restrict pme) 
{
  // unpack parameters
  const real*   zs = (real*) spline->zs;
  const real*   zx = (real*) spline->theta[0];
  const real*   zy = (real*) spline->theta[1];
  const real*   zz = (real*) spline->theta[2];
  const real*   zfx= (real*) spline->dtheta[0];
  const real*   zfy= (real*) spline->dtheta[1];
  const real*   zfz= (real*) spline->dtheta[2];

  const int    p = params->P;
  const double h = params->h;

  int i, j, k, m, mm, idx_zs, idx_zz;
  real force_m[3], Hzc, qm, h3=h*h*h;
  int * idxptr;
  int i0,j0,k0;
  int index_x, index_xy;

#ifdef CALC_ENERGY
  real phi_m;
#endif

  int pny   = pme->pmegrid_ny;
  int pnz   = pme->pmegrid_nz;

/*
#ifdef _OPENMP
#pragma omp for private(mm) schedule(static) // work-share over OpenMP threads here
#endif
*/
  for(mm=0; mm<spline->n; mm++)
    {
      m   = spline->ind[mm];
      qm  = q[m];
      idxptr = atc->idx[mm];
      i0 = idxptr[XX];
      j0 = idxptr[YY];
      k0 = idxptr[ZZ];

      force_m[0] = 0; force_m[1] = 0; force_m[2] = 0;
#ifdef CALC_ENERGY
      phi_m = 0;
#endif

      idx_zs = 0;
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	real zxi = zx[m*p+i];
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  real zxzy = zxi*zy[m*p+j];
	  idx_zz=p*m;
	  for(k = 0; k<p; k++)
	    {
	      Hzc    = grid[index_xy + k]*zs[idx_zs]*zz[idx_zz]*zxzy;
#ifdef CALC_ENERGY
	      phi_m += Hzc;
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

      float factor = -qm*scale*h3;
      force[m][XX] += factor*force_m[0];
      force[m][YY] += factor*force_m[1];
      force[m][ZZ] += factor*force_m[2];

#ifdef CALC_ENERGY
      st->phi[m]   = -h3*scale*phi_m;
#endif
    }
}


// -----------------------------------------------------------------------------
static void 
SE_int_split_SSE_d(rvec * gmx_restrict force,  real * gmx_restrict grid, 
		   real * gmx_restrict q,
		   splinedata_t         * gmx_restrict spline,
		   const SE_params  * gmx_restrict params, 
		   real scale, gmx_bool bClearF,
		   const pme_atomcomm_t * gmx_restrict atc,
		   const gmx_pme_t      * gmx_restrict pme)
{
  // unpack params
  const double*   zs = (double*) spline->zs;
  const double*   zx = (double*) spline->theta[0];
  const double*   zy = (double*) spline->theta[1];
  const double*   zz = (double*) spline->theta[2];
  const double*   zfx= (double*) spline->dtheta[0];
  const double*   zfy= (double*) spline->dtheta[1];
  const double*   zfz= (double*) spline->dtheta[2];


  const int  p = params->P;
  const int  N = params->N;
  const double h = params->h;

  int i,j,k,m,idx_zs,idx_zz, mm,kmp;
  double qm, h3 = h*h*h;
  int * idxptr;
  int i0,j0,k0;
  int index_x, index_xy;
  
  int pny   = pme->pmegrid_ny;
  int pnz   = pme->pmegrid_nz;

  __m128d rH0, rZZ0, rZS0,rZFZ0;
  __m128d rC, rCX, rCY;
  __m128d rFX, rFY, rFZ;
#ifdef CALC_ENERGY
  __m128d rP, rCP;
#endif

  for(mm=0; mm<N; mm++){
    m = spline->ind[mm];
    qm = q[m];
    idxptr = atc->idx[mm];
    i0 = idxptr[XX];
    j0 = idxptr[YY];
    k0 = idxptr[ZZ];
    
    idx_zs = 0;
    rFX = _mm_setzero_pd();
    rFY = _mm_setzero_pd();
    rFZ = _mm_setzero_pd();
    
#ifdef CALC_ENERGY
    rP = _mm_setzero_pd();
#endif
    if(0){ // H[idx] is 16-aligned so vectorization simple
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	real zxi = zx[m*p+i];
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  real zxzy = zxi*zy[m*p+j];

	  rC  = _mm_set1_pd( zxzy );
	  rCX = _mm_set1_pd( zxzy * zfx[m*p + i] );
	  rCY = _mm_set1_pd( zxzy * zfy[m*p + j] );
#ifdef CALC_ENERGY
	  rCP  = _mm_set1_pd( zxzy );
#endif  
	  idx_zz=m*p;
	  for(k = 0; k<p; k+=2){
	    kmp = m*p+k;
	    rZFZ0= _mm_load_pd( zfz+ kmp );
	    rH0  = _mm_load_pd( grid + index_xy + k);
	    rZZ0 = _mm_load_pd( zz + idx_zz);
	    rZS0 = _mm_load_pd( zs + idx_zs);
	    
	    rH0 = _mm_mul_pd(rH0,_mm_mul_pd(rZZ0,rZS0));
	    rFX = _mm_add_pd(rFX,_mm_mul_pd(rH0,rCX));
	    rFY = _mm_add_pd(rFY,_mm_mul_pd(rH0,rCY));
	    rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(rC,rZFZ0)));
	    
#ifdef CALC_ENERGY
	    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,rCP));
#endif			
	    idx_zs+=2; 
	    idx_zz+=2;
	  }
	}
      }
    }
    else { // H[idx] not 16-aligned, so use non-aligned loads
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	real zxi = zx[m*p+i];
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  real zxzy = zxi*zy[m*p+j];

	  rC  = _mm_set1_pd( zxzy );
	  rCX = _mm_set1_pd( zxzy * zfx[m*p + i] );
	  rCY = _mm_set1_pd( zxzy * zfy[m*p + j] );
#ifdef CALC_ENERGY
	  rCP  = _mm_set1_pd( zxzy );
#endif  
	  idx_zz=m*p;
	  for(k = 0; k<p; k+=2){
	    kmp = m*p+k;
	    rZFZ0= _mm_load_pd( zfz+ kmp );
	    rH0  = _mm_loadu_pd( grid + index_xy + k);
	    rZZ0 = _mm_load_pd( zz + idx_zz);
	    rZS0 = _mm_load_pd( zs + idx_zs);
	    
	    rH0 = _mm_mul_pd(rH0,_mm_mul_pd(rZZ0,rZS0));
	    rFX = _mm_add_pd(rFX,_mm_mul_pd(rH0,rCX));
	    rFY = _mm_add_pd(rFY,_mm_mul_pd(rH0,rCY));
	    rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(rC,rZFZ0)));
	    
#ifdef CALC_ENERGY
	    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,rCP));
#endif			
	    idx_zs+=2; 
	    idx_zz+=2;
	  }
	}
      }
    }

    if(bClearF){
      force[m][XX] = 0;
      force[m][YY] = 0;
      force[m][ZZ] = 0;
    }
    
    double factor = -qm*scale*h3;
    force[m][XX] += factor*reduce_sse(rFX);
    force[m][YY] += factor*reduce_sse(rFY);
    force[m][ZZ] += factor*reduce_sse(rFZ);
    
#ifdef CALC_ENERGY
    st->phi[m] = -scale*h3*reduce_sse(rP);
#endif
  }
}


// -----------------------------------------------------------------------------
static void 
SE_int_split_SSE_u8_d(rvec * gmx_restrict force,  real * gmx_restrict grid, 
		      real * gmx_restrict q,
		      splinedata_t         * gmx_restrict spline,
		      const SE_params  * gmx_restrict params, 
		      real scale, gmx_bool bClearF,
		      const pme_atomcomm_t * gmx_restrict atc,
		      const gmx_pme_t      * gmx_restrict pme)
{
  // unpack params
  const double*   zs = (double*) spline->zs;
  const double*   zx = (double*) spline->theta[0];
  const double*   zy = (double*) spline->theta[1];
  const double*   zz = (double*) spline->theta[2];
  const double*   zfx= (double*) spline->dtheta[0];
  const double*   zfy= (double*) spline->dtheta[1];
  const double*   zfz= (double*) spline->dtheta[2];

  const int  p = params->P;
  const int  N = params->N;
  const double h = params->h;

  int i,j,k,m,idx,idx_zs,idx_zz, mm;
  double qm, h3=h*h*h;
  int * idxptr;
  int i0,j0,k0;
  int index_x, index_xy;

  int pny   = pme->pmegrid_ny;
  int pnz   = pme->pmegrid_nz;

  __m128d rH0, rZZ0, rZS0, rZFZ0;
  __m128d rH1, rZZ1, rZS1, rZFZ1;
  __m128d rH2, rZZ2, rZS2, rZFZ2;
  __m128d rH3, rZZ3, rZS3, rZFZ3;
  __m128d rFX, rFY, rFZ;
  __m128d  rC, rCX, rCY;

#ifdef CALC_ENERGY
  __m128d rP, rCP;
#endif

  for(m=0; m<N; m++){
    mm  = spline->ind[m];
    qm = q[mm];
    idxptr = atc->idx[mm];
    i0 = idxptr[XX];
    j0 = idxptr[YY];
    k0 = idxptr[ZZ];
    idx = i0*pny*pnz+j0*pnz+k0;
    
    _mm_prefetch( (void*) (grid+idx), _MM_HINT_T0);
    
    idx_zs = 0;
    _mm_prefetch( (void*) zs, _MM_HINT_T0);
    
    rFX = _mm_setzero_pd();
    rFY = _mm_setzero_pd();
    rFZ = _mm_setzero_pd();
#ifdef CALC_ENERGY
    rP  = _mm_setzero_pd();
#endif
    
    if(0){ // H[idx] is 16-aligned so vectorization simple
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	real zxi = zx[m*p+i];
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  real zxzy = zxi*zy[m*p+j];
	  rC  = _mm_set1_pd( zxzy );
	  rCX = _mm_set1_pd(zxzy*zfx[m*p+i]);
	  rCY = _mm_set1_pd(zxzy*zfy[m*p+j]);
#ifdef CALC_ENERGY
	  rCP = _mm_set1_pd( zxzy );
#endif
	  
	  idx_zz=m*p;
	  
	  for(k = 0; k<p; k+=8){
	    rH0  = _mm_load_pd( grid + index_xy + k    );
	    rH1  = _mm_load_pd( grid + index_xy + k + 2);
	    rH2  = _mm_load_pd( grid + index_xy + k + 4);
	    rH3  = _mm_load_pd( grid + index_xy + k + 6);
	    
	    rZZ0 = _mm_load_pd( zz + idx_zz    );
	    rZZ1 = _mm_load_pd( zz + idx_zz + 2);
	    rZZ2 = _mm_load_pd( zz + idx_zz + 4);
	    rZZ3 = _mm_load_pd( zz + idx_zz + 6);
	    
	    rZS0 = _mm_load_pd( zs + idx_zs    );
	    rZS1 = _mm_load_pd( zs + idx_zs + 2);
	    rZS2 = _mm_load_pd( zs + idx_zs + 4);
	    rZS3 = _mm_load_pd( zs + idx_zs + 6);
	    
	    rZFZ0 = _mm_load_pd(zfz+ idx_zz    );
	    rZFZ1 = _mm_load_pd(zfz+ idx_zz + 2);
	    rZFZ2 = _mm_load_pd(zfz+ idx_zz + 4);
	    rZFZ3 = _mm_load_pd(zfz+ idx_zz + 6);

	    __m128d fac0 = _mm_mul_pd(rH0, _mm_mul_pd(rZZ0, rZS0));
	    __m128d fac1 = _mm_mul_pd(rH1, _mm_mul_pd(rZZ1, rZS1));
	    __m128d fac2 = _mm_mul_pd(rH2, _mm_mul_pd(rZZ2, rZS2));
	    __m128d fac3 = _mm_mul_pd(rH3, _mm_mul_pd(rZZ3, rZS3));
	    
	    __m128d fac  = _mm_add_pd(_mm_add_pd(fac0,fac1),_mm_add_pd(fac2,fac3));

	    rFX = _mm_add_pd(rFX,_mm_mul_pd(fac,rCX));
	    rFY = _mm_add_pd(rFY,_mm_mul_pd(fac,rCY));
	    
	    fac0 = _mm_mul_pd(fac0,rZFZ0);
	    fac1 = _mm_mul_pd(fac1,rZFZ1);
	    fac2 = _mm_mul_pd(fac2,rZFZ2);
	    fac3 = _mm_mul_pd(fac3,rZFZ3);

	    rFZ  = _mm_add_pd(rFZ,_mm_mul_pd(_mm_add_pd(fac0,fac1),rC));

#ifdef CALC_ENERGY
	    rP =_mm_add_pd(rP,_mm_mul_pd(fac,rC));
#endif
	    idx_zs+=8; 
	    idx_zz+=8;
	  }
	}
      }
    }
    else{ // H[idx] not 16-aligned, so use non-aligned load from H
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	real zxi = zx[m*p+i];
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  real zxzy = zxi*zy[m*p+j];
	  rC  = _mm_set1_pd( zxzy );
	  rCX = _mm_set1_pd(zxzy*zfx[m*p+i]);
	  rCY = _mm_set1_pd(zxzy*zfy[m*p+j]);
#ifdef CALC_ENERGY
	  rCP = _mm_set1_pd( zxzy );
#endif
	  
	  idx_zz=m*p;
	  
	  for(k = 0; k<p; k+=8){
	    rH0  = _mm_loadu_pd( grid + index_xy + k    );
	    rH1  = _mm_loadu_pd( grid + index_xy + k + 2);
	    rH2  = _mm_loadu_pd( grid + index_xy + k + 4);
	    rH3  = _mm_loadu_pd( grid + index_xy + k + 6);
	    
	    rZZ0 = _mm_load_pd( zz + idx_zz    );
	    rZZ1 = _mm_load_pd( zz + idx_zz + 2);
	    rZZ2 = _mm_load_pd( zz + idx_zz + 4);
	    rZZ3 = _mm_load_pd( zz + idx_zz + 6);
	    
	    rZS0 = _mm_load_pd( zs + idx_zs    );
	    rZS1 = _mm_load_pd( zs + idx_zs + 2);
	    rZS2 = _mm_load_pd( zs + idx_zs + 4);
	    rZS3 = _mm_load_pd( zs + idx_zs + 6);
	    
	    rZFZ0 = _mm_load_pd(zfz+ idx_zz    );
	    rZFZ1 = _mm_load_pd(zfz+ idx_zz + 2);
	    rZFZ2 = _mm_load_pd(zfz+ idx_zz + 4);
	    rZFZ3 = _mm_load_pd(zfz+ idx_zz + 6);

	    __m128d fac0 = _mm_mul_pd(rH0, _mm_mul_pd(rZZ0, rZS0));
	    __m128d fac1 = _mm_mul_pd(rH1, _mm_mul_pd(rZZ1, rZS1));
	    __m128d fac2 = _mm_mul_pd(rH2, _mm_mul_pd(rZZ2, rZS2));
	    __m128d fac3 = _mm_mul_pd(rH3, _mm_mul_pd(rZZ3, rZS3));
	    
	    __m128d fac  = _mm_add_pd(_mm_add_pd(fac0,fac1),_mm_add_pd(fac2,fac3));

	    rFX = _mm_add_pd(rFX,_mm_mul_pd(fac,rCX));
	    rFY = _mm_add_pd(rFY,_mm_mul_pd(fac,rCY));
	    
	    fac0 = _mm_mul_pd(fac0,rZFZ0);
	    fac1 = _mm_mul_pd(fac1,rZFZ1);
	    fac2 = _mm_mul_pd(fac2,rZFZ2);
	    fac3 = _mm_mul_pd(fac3,rZFZ3);

	    rFZ  = _mm_add_pd(rFZ,_mm_mul_pd(_mm_add_pd(fac0,fac1),rC));

#ifdef CALC_ENERGY
	    rP =_mm_add_pd(rP,_mm_mul_pd(fac,rC));
#endif
	    idx_zs+=8; 
	    idx_zz+=8;
	  }
	}
      }
    }

    if(bClearF){
      force[m][XX] = 0;
      force[m][YY] = 0;
      force[m][ZZ] = 0;
    }
    
    double factor = -qm*scale*h3;
    force[m][XX] += factor*reduce_sse(rFX);
    force[m][YY] += factor*reduce_sse(rFY);
    force[m][ZZ] += factor*reduce_sse(rFZ);
    
#ifdef CALC_ENERGY
    _mm_store_pd(s,rP);
    st->phi[m] = -scale*h3*reduce_sse(rP);
#endif
  }
}

// -----------------------------------------------------------------------------
static void 
SE_int_split_SSE_P8_d(rvec * gmx_restrict force, real * gmx_restrict grid, 
		      real * gmx_restrict q,
		      splinedata_t          * gmx_restrict spline,
		      const SE_params   * gmx_restrict params, 
		      real scale, gmx_bool bClearF,
		      const pme_atomcomm_t  * gmx_restrict atc,
		      const gmx_pme_t       * gmx_restrict pme)
{
  // unpack params
  const double*   zs = (double*) spline->zs;
  const double*   zx = (double*) spline->theta[0];
  const double*   zy = (double*) spline->theta[1];
  const double*   zz = (double*) spline->theta[2];
  const double*   zfx= (double*) spline->dtheta[0];
  const double*   zfy= (double*) spline->dtheta[1];
  const double*   zfz= (double*) spline->dtheta[2];

  /* ASSUME P=8 const int p = params->P; */
  const int  N = params->N;
  const double h = params->h;

  int i,j, m, idx_zs, mm,m8;
  double qm, h3=h*h*h;
  int * idxptr;
  int i0,j0,k0;
  int index_x, index_xy;

  int pny   = pme->pmegrid_ny;
  int pnz   = pme->pmegrid_nz;

  // hold entire zz vector
  __m128d rZZ0, rZZ1, rZZ2, rZZ3;
  __m128d rC, rCX, rCY;
  __m128d rH0, rH1, rH2, rH3; 
  __m128d rZS0, rZS1, rZS2, rZS3;
  __m128d rZFZ0, rZFZ1, rZFZ2, rZFZ3;
  __m128d rFX, rFY, rFZ;

#ifdef CALC_ENERGY
  __m128d rP, rCP;
#endif

  for(mm=0; mm<N; mm++){
    m = spline->ind[mm];
    qm = q[m];
    idxptr = atc->idx[mm];
    i0 = idxptr[XX];
    j0 = idxptr[YY];
    k0 = idxptr[ZZ];
    m8 = m*8;
    
    idx_zs = 0;
    rFX = _mm_setzero_pd();
    rFY = _mm_setzero_pd();
    rFZ = _mm_setzero_pd();
#ifdef CALC_ENERGY
    rP  = _mm_setzero_pd();
#endif
    
    /* hoist load of ZZ vector */
    rZZ0 = _mm_load_pd(zz + m8     );
    rZZ1 = _mm_load_pd(zz + m8 + 2 );
    rZZ2 = _mm_load_pd(zz + m8 + 4 );
    rZZ3 = _mm_load_pd(zz + m8 + 6 );
    
    /* hoist load of ZFZ vector */
    rZFZ0 = _mm_load_pd(zfz + m8     );
    rZFZ1 = _mm_load_pd(zfz + m8 + 2 );
    rZFZ2 = _mm_load_pd(zfz + m8 + 4 );
    rZFZ3 = _mm_load_pd(zfz + m8 + 6 );
    
    if(0){ // H[idx] is 16-aligned so vectorization simple
      for(i = 0; i<8; i++){
	index_x = (i0+i)*pny*pnz;
	real zxi = zx[m8+i];
	for(j = 0; j<8; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  real zxzy = zxi*zy[m8+j];
	  rC  = _mm_set1_pd( zxzy );
	  rCX = _mm_set1_pd( zxzy * zfx[m8+i]);
	  rCY = _mm_set1_pd( zxzy * zfy[m8+j]);
#ifdef CALC_ENERGY
	  rCP  = _mm_set1_pd( zxzy );
#endif
	  
	  rH0  = _mm_load_pd( grid + index_xy    );
	  rH1  = _mm_load_pd( grid + index_xy + 2);
	  rH2  = _mm_load_pd( grid + index_xy + 4);
	  rH3  = _mm_load_pd( grid + index_xy + 6);
		  
	  rZS0 = _mm_load_pd( zs + idx_zs    );
	  rZS1 = _mm_load_pd( zs + idx_zs + 2);
	  rZS2 = _mm_load_pd( zs + idx_zs + 4);
	  rZS3 = _mm_load_pd( zs + idx_zs + 6);

	  rH0 = _mm_mul_pd(rH0,_mm_mul_pd(rZZ0,rZS0));
	  rH1 = _mm_mul_pd(rH1,_mm_mul_pd(rZZ1,rZS1));
	  rH2 = _mm_mul_pd(rH2,_mm_mul_pd(rZZ2,rZS2));
	  rH3 = _mm_mul_pd(rH3,_mm_mul_pd(rZZ3,rZS3));

	  __m128d rh = _mm_add_pd(_mm_add_pd(rH0,rH1),_mm_add_pd(rH2,rH3));

	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rh,rCX));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rh,rCY));
	  
	  __m128d fac0 = _mm_mul_pd(rH0,rZFZ0);
	  __m128d fac1 = _mm_mul_pd(rH1,rZFZ1);
	  __m128d fac2 = _mm_mul_pd(rH2,rZFZ2);
	  __m128d fac3 = _mm_mul_pd(rH3,rZFZ3);
	  rh = _mm_add_pd(_mm_add_pd(fac0,fac1),_mm_add_pd(fac2,fac3));

	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rh,rC));

#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rh,rCP));
#endif

	  idx_zs +=8;
	}
      }
    }
    else{ // H[idx] not 16-aligned, so use non-aligned loads
      for(i = 0; i<8; i++){
	index_x = (i0+i)*pny*pnz;
	real zxi = zx[m8+i];
	for(j = 0; j<8; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  real zxzy = zxi*zy[m8+j];
	  rC  = _mm_set1_pd( zxzy );
	  rCX = _mm_set1_pd( zxzy * zfx[m8+i]);
	  rCY = _mm_set1_pd( zxzy * zfy[m8+j]);
#ifdef CALC_ENERGY
	  rCP  = _mm_set1_pd( zxzy );
#endif
	  
	  rH0  = _mm_loadu_pd( grid + index_xy    );
	  rH1  = _mm_loadu_pd( grid + index_xy + 2);
	  rH2  = _mm_loadu_pd( grid + index_xy + 4);
	  rH3  = _mm_loadu_pd( grid + index_xy + 6);
		  
	  rZS0 = _mm_load_pd( zs + idx_zs    );
	  rZS1 = _mm_load_pd( zs + idx_zs + 2);
	  rZS2 = _mm_load_pd( zs + idx_zs + 4);
	  rZS3 = _mm_load_pd( zs + idx_zs + 6);

	  rH0 = _mm_mul_pd(rH0,_mm_mul_pd(rZZ0,rZS0));
	  rH1 = _mm_mul_pd(rH1,_mm_mul_pd(rZZ1,rZS1));
	  rH2 = _mm_mul_pd(rH2,_mm_mul_pd(rZZ2,rZS2));
	  rH3 = _mm_mul_pd(rH3,_mm_mul_pd(rZZ3,rZS3));

	  __m128d rh = _mm_add_pd(_mm_add_pd(rH0,rH1),_mm_add_pd(rH2,rH3));

	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rh,rCX));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rh,rCY));
	  
	  __m128d fac0 = _mm_mul_pd(rH0,rZFZ0);
	  __m128d fac1 = _mm_mul_pd(rH1,rZFZ1);
	  __m128d fac2 = _mm_mul_pd(rH2,rZFZ2);
	  __m128d fac3 = _mm_mul_pd(rH3,rZFZ3);
	  rh = _mm_add_pd(_mm_add_pd(fac0,fac1),_mm_add_pd(fac2,fac3));

	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rh,rC));

#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rh,rCP));
#endif

	  idx_zs +=8;
	}
      }      
    }

    if(bClearF){
      force[m][XX] = 0;
      force[m][YY] = 0;
      force[m][ZZ] = 0;
    }

    double factor = -qm*scale*h3;
    force[m][XX] += factor*reduce_sse(rFX);
    force[m][YY] += factor*reduce_sse(rFY);
    force[m][ZZ] += factor*reduce_sse(rFZ);
    
#ifdef CALC_ENERGY
    _mm_store_pd(s,rP);
    st->phi[m] = -scale*h3*reduce_sse(rP);
#endif
  }
}

// -----------------------------------------------------------------------------
static void 
SE_int_split_SSE_P16_d(rvec * gmx_restrict force, real * gmx_restrict grid, 
		       real * gmx_restrict q,
		       splinedata_t          * gmx_restrict spline,
		       const SE_params   * gmx_restrict params, 
		       real scale, gmx_bool bClearF,
		       const pme_atomcomm_t  * gmx_restrict atc,
		       const gmx_pme_t       * gmx_restrict pme)
{
  // unpack params
  const double*   zs = (double*) spline->zs;
  const double*   zx = (double*) spline->theta[0];
  const double*   zy = (double*) spline->theta[1];
  const double*   zz = (double*) spline->theta[2];
  const double*   zfx= (double*) spline->dtheta[0];
  const double*   zfy= (double*) spline->dtheta[1];
  const double*   zfz= (double*) spline->dtheta[2];

  /* ASSUME P=16 const int p = params->P; */
  const int  N = params->N;
  const double h = params->h;

  int i,j,m,idx,idx_zs, mm, m16,im16,jm16;
  double qm, h3=h*h*h;
  int * idxptr;
  int i0,j0,k0;
  int index_x, index_xy;

  int pny   = pme->pmegrid_ny;
  int pnz   = pme->pmegrid_nz;

  // hold entire zz vector
  __m128d rZZ0 , rZZ1 , rZZ2 , rZZ3 , rZZ4 , rZZ5 , rZZ6 , rZZ7; 
  __m128d rZFZ0, rZFZ1, rZFZ2, rZFZ3, rZFZ4, rZFZ5, rZFZ6, rZFZ7;
  __m128d rC, rCX, rCY, rFX, rFY, rFZ;
  __m128d rH0, rZS0;

#ifdef CALC_ENERGY
  __m128d rP, rCP;
#endif

  for(mm=0; mm<N; mm++){
    m = spline->ind[mm];
    qm = q[m];
    idxptr = atc->idx[mm];
    i0 = idxptr[XX];
    j0 = idxptr[YY];
    k0 = idxptr[ZZ];
    idx = i0*pny*pnz+j0*pnz+k0;
    m16 = m*16;
    
    _mm_prefetch( (void*) (grid+idx), _MM_HINT_T0);
    
    idx_zs = 0;
    _mm_prefetch( (void*) zs, _MM_HINT_T0);
    
    rFX = _mm_setzero_pd();
    rFY = _mm_setzero_pd();
    rFZ = _mm_setzero_pd();
#ifdef CALC_ENERGY
    rP  = _mm_setzero_pd();
#endif
    
    /* hoist load of ZZ vector */
    rZZ0 = _mm_load_pd(zz + m16     );
    rZZ1 = _mm_load_pd(zz + m16 + 2 );
    rZZ2 = _mm_load_pd(zz + m16 + 4 );
    rZZ3 = _mm_load_pd(zz + m16 + 6 );
    rZZ4 = _mm_load_pd(zz + m16 + 8 );
    rZZ5 = _mm_load_pd(zz + m16 + 10);
    rZZ6 = _mm_load_pd(zz + m16 + 12);
    rZZ7 = _mm_load_pd(zz + m16 + 14);

    /* hoist load of ZFZ vector */
    rZFZ0 = _mm_load_pd(zfz + m16     );
    rZFZ1 = _mm_load_pd(zfz + m16 + 2 );
    rZFZ2 = _mm_load_pd(zfz + m16 + 4 );
    rZFZ3 = _mm_load_pd(zfz + m16 + 6 );
    rZFZ4 = _mm_load_pd(zfz + m16 + 8 );
    rZFZ5 = _mm_load_pd(zfz + m16 + 10);
    rZFZ6 = _mm_load_pd(zfz + m16 + 12);
    rZFZ7 = _mm_load_pd(zfz + m16 + 14);

    if(0){ // H[idx] is 16-aligned so vectorization simple
      for(i = 0; i<16; i++){
	index_x = (i0+i)*pny*pnz;
	real zxi = zx[m16+i];
	im16 = i + m16;
	for(j = 0; j<16; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  real zxzy = zxi*zy[m16+j];
	  jm16 = j + m16;
	  rC  = _mm_set1_pd( zxzy);
	  rCX = _mm_set1_pd( zxzy*zfx[im16]);
	  rCY = _mm_set1_pd( zxzy*zfy[jm16]);
#ifdef CALC_ENERGY
	  rCP = _mm_set1_pd( zxzy);
#endif
	  
	  /* 0 */ 
	  rH0  = _mm_load_pd( grid + index_xy );
	  rZS0 = _mm_load_pd( zs + idx_zs);		    
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
	  
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ0,rZZ0),rCP),rZS0)));
#endif
	  /* 1 */ 
	  rH0  = _mm_load_pd( grid + index_xy + 2);
	  rZS0 = _mm_load_pd( zs + idx_zs + 2);		    
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ1,rZZ1),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rCP),rZS0)));
#endif
	  
	  /* 2 */ 
	  rH0  = _mm_load_pd( grid + index_xy + 4);
	  rZS0 = _mm_load_pd( zs + idx_zs + 4);		    
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ2,rZZ2),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rCP),rZS0)));
#endif
	  
	  /* 3 */ 
	  rH0  = _mm_load_pd( grid + index_xy + 6);
	  rZS0 = _mm_load_pd( zs + idx_zs + 6);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ3,rZZ3),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rCP),rZS0)));
#endif

	  /* 4 */ 
	  rH0  = _mm_load_pd( grid + index_xy + 8);
	  rZS0 = _mm_load_pd( zs + idx_zs + 8);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ4,rZZ4),rC),rZS0)));
#ifdef CALC_NERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rCP),rZS0)));
#endif

	  /* 5 */ 
	  rH0  = _mm_load_pd( grid + index_xy + 10);
	  rZS0 = _mm_load_pd( zs + idx_zs + 10);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ5,rZZ5),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rCP),rZS0)));
#endif

	  /* 6 */ 
	  rH0  = _mm_load_pd( grid + index_xy + 12);
	  rZS0 = _mm_load_pd( zs + idx_zs + 12);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ6,rZZ6),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rCP),rZS0)));
#endif
	  /* 7 */ 
	  rH0  = _mm_load_pd( grid + index_xy + 14);
	  rZS0 = _mm_load_pd( zs + idx_zs + 14);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ7,rZZ7),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rCP),rZS0)));
#endif

	  idx_zs +=16;
	}
      }
    }
    else{ // H[idx] not 16-aligned, so use non-aligned loads
      for(i = 0; i<16; i++){
	index_x = (i0+i)*pny*pnz;
	real zxi = zx[m16+i];
	im16 = i + m16;
	for(j = 0; j<16; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  real zxzy = zxi*zy[m16+j];
	  jm16 = j + m16;
	  rC  = _mm_set1_pd( zxzy);
	  rCX = _mm_set1_pd( zxzy*zfx[im16]);
	  rCY = _mm_set1_pd( zxzy*zfy[jm16]);
#ifdef CALC_ENERGY
	  rCP = _mm_set1_pd( zxzy);
#endif
	  
	  /* 0 */ 
	  rH0  = _mm_load_pd( grid + index_xy );
	  rZS0 = _mm_load_pd( zs + idx_zs);		    
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
	  
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ0,rZZ0),rCP),rZS0)));
#endif
	  /* 1 */ 
	  rH0  = _mm_load_pd( grid + index_xy + 2);
	  rZS0 = _mm_load_pd( zs + idx_zs + 2);		    
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ1,rZZ1),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rCP),rZS0)));
#endif
	  
	  /* 2 */ 
	  rH0  = _mm_load_pd( grid + index_xy + 4);
	  rZS0 = _mm_load_pd( zs + idx_zs + 4);		    
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ2,rZZ2),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rCP),rZS0)));
#endif
	  
	  /* 3 */ 
	  rH0  = _mm_load_pd( grid + index_xy + 6);
	  rZS0 = _mm_load_pd( zs + idx_zs + 6);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ3,rZZ3),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rCP),rZS0)));
#endif

	  /* 4 */ 
	  rH0  = _mm_load_pd( grid + index_xy + 8);
	  rZS0 = _mm_load_pd( zs + idx_zs + 8);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ4,rZZ4),rC),rZS0)));
#ifdef CALC_NERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rCP),rZS0)));
#endif

	  /* 5 */ 
	  rH0  = _mm_load_pd( grid + index_xy + 10);
	  rZS0 = _mm_load_pd( zs + idx_zs + 10);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ5,rZZ5),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rCP),rZS0)));
#endif

	  /* 6 */ 
	  rH0  = _mm_load_pd( grid + index_xy + 12);
	  rZS0 = _mm_load_pd( zs + idx_zs + 12);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ6,rZZ6),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rCP),rZS0)));
#endif
	  /* 7 */ 
	  rH0  = _mm_load_pd( grid + index_xy + 14);
	  rZS0 = _mm_load_pd( zs + idx_zs + 14);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ7,rZZ7),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rCP),rZS0)));
#endif

	  idx_zs +=16;
	}
      }      
    }

    if(bClearF){
      force[m][XX] = 0;
      force[m][YY] = 0;
      force[m][ZZ] = 0;
    }

    double factor = -qm*scale*h3;
    force[m][XX] += factor*reduce_sse(rFX);
    force[m][YY] += factor*reduce_sse(rFY);
    force[m][ZZ] += factor*reduce_sse(rFZ);
    
#ifdef CALC_ENERGY
    _mm_store_pd(s,rP);
    st->phi[m] = -scale*h3*reduce_sse(rP);
#endif
  }
}

// -----------------------------------------------------------------------------
static void 
SE_int_split_SSE_dispatch_d(rvec* force, real* grid, real* q,
			  splinedata_t *spline,
			  const SE_params* params, real scale,
			  gmx_bool bClearF,
			  const pme_atomcomm_t *atc,
			  const gmx_pme_t *pme)
{
  const int p = params->P;
  const int incrj = params->dims[2]; // middle increment
  const int incri = params->npdims[2]*(params->dims[1]);// outer increment

  // if P is odd, or if either increment is odd, fall back on vanilla
  if( is_odd(p) || is_odd(incri) || is_odd(incrj) ){
    __DISPATCHER_MSG("[FGG INT SSE DOUBLE] SSE Abort (PARAMS)\n");
    SE_int_split_d(force, grid, q, spline, params, scale, bClearF, atc, pme);
    return;
  }
  // otherwise the preconditions for SSE codes are satisfied. 
  if(p==8){
    // specific for p=8
    __DISPATCHER_MSG("[FGG INT SSE DOUBLE] P=8\n");
    SE_int_split_SSE_P8_d(force, grid, q, spline, params, scale, bClearF, atc, pme);
  }
  else if(p==16){
    // specific for p=16
    __DISPATCHER_MSG("[FGG INT SSE DOUBLE] P=16\n");
    SE_int_split_SSE_P16_d(force, grid, q, spline, params, scale, bClearF, atc, pme); 
  }
  else if(p%8==0){
    // for p divisible by 8
    __DISPATCHER_MSG("[FGG INT SSE DOUBLE] P unroll 8\n");
    SE_int_split_SSE_u8_d(force, grid, q, spline, params, scale, bClearF, atc, pme); 
  }
  else{
    // vanilla SSE code (any even p)
    __DISPATCHER_MSG("[FGG INT SSE DOUBLE] Vanilla\n");
    SE_int_split_SSE_d(force, grid, q, spline, params, scale, bClearF, atc, pme);
  }
}


#endif // GMX_DOUBLE
#endif // _SE_INT_SSE_DOUBLE_
