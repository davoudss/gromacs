#ifndef _SE_INT_AVX_256_DOUBLE_
#define _SE_INT_AVX_256_DOUBLE_


/* SE AVX 256 double integration */
#if GMX_DOUBLE==1
#if GMX_SIMD_X86_AVX_256
#include "se.h"

static inline double
reduce(__m256d a)
{
  __m128d V = _mm256_castpd256_pd128(a);
  __m128d W = _mm256_extractf128_pd(a, 0x1);
  V = _mm_add_pd(W,V);
  V = _mm_hadd_pd(V,V);
  return *reinterpret_cast<double *>(&V);
}

static void
SE_int_split_AVX_gaussian_d(rvec * gmx_restrict force, real * gmx_restrict grid, 
			    real * gmx_restrict q,
			    splinedata_t          * gmx_restrict spline,
			    const SE_params   * gmx_restrict params, 
			    real scale, gmx_bool bClearF,
			    const pme_atomcomm_t  * gmx_restrict atc,
			    const gmx_pme_t       * gmx_restrict pme)
{
  // unpack params
  const double* zs = (double*) spline->zs;
  const double* zx = (double*) spline->theta[0];
  const double* zy = (double*) spline->theta[1];
  const double* zz = (double*) spline->theta[2];
  const double* zfx =(double*) spline->dtheta[0];
  const double* zfy =(double*) spline->dtheta[1];
  const double* zfz =(double*) spline->dtheta[2];

  const int    p = params->P;
  const int    N = params->N;
  const double h = params->h;

  int i,j,k,m, idx_zs,idx_zz, mm;
  double qm, h3=h*h*h;
  int * idxptr;
  int i0,j0,k0;
  int index_x, index_xy;

  int pny   = pme->pmegrid_ny;
  int pnz   = pme->pmegrid_nz;

  __m256d rH0, rZZ0, rZS0,rZFZ0;
  __m256d rC, rCX, rCY;
  __m256d rFX, rFY, rFZ;
#ifdef CALC_ENERGY
  __m256d rP, rCP;
#endif

  for(m=0; m<N; m++){
    mm  = spline->ind[m];
    qm  = q[mm];
    idxptr = atc->idx[mm];
    i0 = idxptr[XX];
    j0 = idxptr[YY];
    k0 = idxptr[ZZ];

    idx_zs = 0;
    rFX = _mm256_setzero_pd();
    rFY = _mm256_setzero_pd();
    rFZ = _mm256_setzero_pd();
#ifdef CALC_ENERGY
    rP = _mm256_setzero_pd();
#endif
    
    if(0){ // H[idx] is 32-aligned so vectorization simple
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	real zxi = zx[m*p+i];
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  real zxzy = zxi*zy[m*p+j];
	  rC  = _mm256_set1_pd( zxzy );
	  rCX = _mm256_set1_pd( zxzy*zfx[m*p+i] );
	  rCY = _mm256_set1_pd( zxzy*zfy[m*p+j] );
#ifdef CALC_ENERGY
	  rCP  = _mm256_set1_pd( zxzy );
#endif

	  idx_zz=m*p;
	  for(k = 0; k<p; k+=4){
	    rZFZ0= _mm256_load_pd( zfz+ m*p+k );
	    rH0  = _mm256_load_pd( grid + index_xy + k );
	    rZZ0 = _mm256_load_pd( zz + idx_zz);
	    rZS0 = _mm256_load_pd( zs + idx_zs);

	    __m256d fac = _mm256_mul_pd(rH0, _mm256_mul_pd(rZZ0, rZS0));
 
	    rFX = _mm256_add_pd(rFX,_mm256_mul_pd(fac,rCX));
	    rFY = _mm256_add_pd(rFY,_mm256_mul_pd(fac,rCY));
	    rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(fac,_mm256_mul_pd(rC,rZFZ0))); 

#ifdef CALC_ENERGY
	    rP =_mm256_add_pd(rP,_mm256_mul_pd(fac,rC));
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
	real zxi = zx[m*p+i];
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  real zxzy = zxi*zy[m*p+j];
	  rC  = _mm256_set1_pd( zxzy );
	  rCX = _mm256_set1_pd( zxzy*zfx[m*p+i] );
	  rCY = _mm256_set1_pd( zxzy*zfy[m*p+j] );
#ifdef CALC_ENERGY
	  rCP  = _mm256_set1_pd( zxzy );
#endif

	  idx_zz=m*p;
	  for(k = 0; k<p; k+=4){
	    rZFZ0= _mm256_load_pd( zfz+ m*p+k );
	    rH0  = _mm256_loadu_pd( grid + index_xy + k);
	    rZZ0 = _mm256_load_pd( zz + idx_zz);
	    rZS0 = _mm256_load_pd( zs + idx_zs);

	    __m256d fac = _mm256_mul_pd(rH0, _mm256_mul_pd(rZZ0, rZS0));
 
	    rFX = _mm256_add_pd(rFX,_mm256_mul_pd(fac,rCX));
	    rFY = _mm256_add_pd(rFY,_mm256_mul_pd(fac,rCY));
	    rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(fac,_mm256_mul_pd(rC,rZFZ0))); 

#ifdef CALC_ENERGY
	    rP =_mm256_add_pd(rP,_mm256_mul_pd(fac,rC));
#endif			
	    idx_zs+=4; 
	    idx_zz+=4;
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
    force[m][XX] += factor*(reduce(rFX));
    force[m][YY] += factor*(reduce(rFY));
    force[m][ZZ] += factor*(reduce(rFZ));

#ifdef CALC_ENERGY
    st->phi[m] = -scale*h3*reduce(rP);
#endif
  }
}

// -----------------------------------------------------------------------------
static void 
SE_int_split_AVX_u8_gaussian_d(rvec * gmx_restrict force, real * gmx_restrict grid, 
			       real * gmx_restrict q,
			       splinedata_t          * gmx_restrict spline,
			       const SE_params   * gmx_restrict params, 
			       real scale, gmx_bool bClearF,
			       const pme_atomcomm_t  * gmx_restrict atc,
			       const gmx_pme_t       * gmx_restrict pme)
{
  // unpack params
  const double*  zs = (double*) spline->zs;
  const double*  zx = (double*) spline->theta[0];
  const double*  zy = (double*) spline->theta[1];
  const double*  zz = (double*) spline->theta[2];
  const double* zfx = (double*) spline->dtheta[0];
  const double* zfy = (double*) spline->dtheta[1];
  const double* zfz = (double*) spline->dtheta[2];

  const int    p = params->P;
  const int    N = params->N;
  const double h = params->h;

  int i,j,k,idx,idx_zs,idx_zz,m, mm;
  double qm, h3=h*h*h;
  int * idxptr;
  int i0,j0,k0;
  int index_x, index_xy;

  int pny   = pme->pmegrid_ny;
  int pnz   = pme->pmegrid_nz;

  __m256d rH0, rZZ0, rZS0, rZFZ0;
  __m256d rH1, rZZ1, rZS1, rZFZ1;
  __m256d rFX, rFY, rFZ;
  __m256d  rC, rCX, rCY;

#ifdef CALC_ENERGY
  __m256d rP, rCP;
#endif

  for(m=0; m<N; m++){
    idx = spline->idx[m];
    mm = spline->ind[m];
    qm = q[mm];
    idxptr = atc->idx[mm];
    i0 = idxptr[XX];
    j0 = idxptr[YY];
    k0 = idxptr[ZZ];
    idx = i0*pny*pnz+j0*pnz+k0;

    _mm_prefetch( (void*) (grid+idx), _MM_HINT_T0);
    idx_zs = 0;
    _mm_prefetch( (void*) zs, _MM_HINT_T0);
    rFX = _mm256_setzero_pd();
    rFY = _mm256_setzero_pd();
    rFZ = _mm256_setzero_pd();
#ifdef CALC_ENERGY
    rP  = _mm256_setzero_pd();
#endif

    if(0){ // H[idx] is 32-aligned so vectorization simple
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	real zxi = zx[m*p+i];
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  real zxzy = zxi*zy[m*p+j];
	  rC  = _mm256_set1_pd( zxzy );
	  rCX = _mm256_set1_pd( zxzy*zfx[m*p+i] );
	  rCY = _mm256_set1_pd( zxzy*zfy[m*p+j] );
#ifdef CALC_ENERGY
	  rCP = _mm256_set1_pd( zxzy );
#endif

	  idx_zz=m*p;
	  for(k = 0; k<p; k+=8){
	    rH0  = _mm256_load_pd( grid + index_xy + k );
	    rH1  = _mm256_load_pd( grid + index_xy + k + 4);

	    rZZ0 = _mm256_load_pd( zz + idx_zz    );
	    rZZ1 = _mm256_load_pd( zz + idx_zz + 4);

	    rZS0 = _mm256_load_pd( zs + idx_zs    );
	    rZS1 = _mm256_load_pd( zs + idx_zs + 4);

	    rZFZ0 = _mm256_load_pd(zfz+ idx_zz    );
	    rZFZ1 = _mm256_load_pd(zfz+ idx_zz + 4);

	    __m256d fac0 = _mm256_mul_pd(rH0, _mm256_mul_pd(rZZ0, rZS0));
	    __m256d fac1 = _mm256_mul_pd(rH1, _mm256_mul_pd(rZZ1, rZS1));
	    __m256d fac = _mm256_add_pd(fac0,fac1);

	    rFX = _mm256_add_pd(rFX,_mm256_mul_pd(fac,rCX));
	    rFY = _mm256_add_pd(rFY,_mm256_mul_pd(fac,rCY));

	    fac0 = _mm256_mul_pd(fac0,rZFZ0);
	    fac1 = _mm256_mul_pd(fac1,rZFZ1);

	    rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(_mm256_add_pd(fac0,fac1),rC));


#ifdef CALC_ENERGY
	    rP =_mm256_add_pd(rP,_mm256_mul_pd(fac,rC));
#endif
	    
	    idx_zs+=8; 
	    idx_zz+=8;
	  }
	}
      }
    }
    else{ // H[idx] not 32-aligned, so use non-aligned load from H
      for(i = 0; i<p; i++){
	index_x = (i0+i)*pny*pnz;
	real zxi = zx[m*p+i];
	for(j = 0; j<p; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  real zxzy = zxi*zy[m*p+j];
	  rC  = _mm256_set1_pd( zxzy );
	  rCX = _mm256_set1_pd( zxzy*zfx[m*p+i] );
	  rCY = _mm256_set1_pd( zxzy*zfy[m*p+j] );
#ifdef CALC_ENERGY
	  rCP = _mm256_set1_pd( zxzy );
#endif

	  idx_zz=m*p;
	  for(k = 0; k<p; k+=8){
	    rH0  = _mm256_loadu_pd( grid + index_xy + k );
	    rH1  = _mm256_loadu_pd( grid + index_xy + k + 4);

	    rZZ0 = _mm256_load_pd( zz + idx_zz    );
	    rZZ1 = _mm256_load_pd( zz + idx_zz + 4);

	    rZS0 = _mm256_load_pd( zs + idx_zs    );
	    rZS1 = _mm256_load_pd( zs + idx_zs + 4);

	    rZFZ0 = _mm256_load_pd(zfz+ idx_zz    );
	    rZFZ1 = _mm256_load_pd(zfz+ idx_zz + 4);

	    __m256d fac0 = _mm256_mul_pd(rH0, _mm256_mul_pd(rZZ0, rZS0));
	    __m256d fac1 = _mm256_mul_pd(rH1, _mm256_mul_pd(rZZ1, rZS1));
	    __m256d fac = _mm256_add_pd(fac0,fac1);

	    rFX = _mm256_add_pd(rFX,_mm256_mul_pd(fac,rCX));
	    rFY = _mm256_add_pd(rFY,_mm256_mul_pd(fac,rCY));

	    fac0 = _mm256_mul_pd(fac0,rZFZ0);
	    fac1 = _mm256_mul_pd(fac1,rZFZ1);

	    rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(_mm256_add_pd(fac0,fac1),rC));

#ifdef CALC_ENERGY
	    rP =_mm256_add_pd(rP,_mm256_mul_pd(fac,rC));
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
    force[m][XX] += factor*(reduce(rFX));
    force[m][YY] += factor*(reduce(rFY));
    force[m][ZZ] += factor*(reduce(rFZ));

#ifdef CALC_ENERGY
    st->phi[m] = -scale*h3*reduce(rP);
#endif
  }
}

// -----------------------------------------------------------------------------
static void
SE_int_split_AVX_P8_gaussian_d(rvec * gmx_restrict force, real * gmx_restrict grid, 
			       real * gmx_restrict q,
			       splinedata_t          * gmx_restrict spline,
			       const SE_params   * gmx_restrict params, 
			       real scale, gmx_bool bClearF,
			       const pme_atomcomm_t  * gmx_restrict atc,
			       const gmx_pme_t       * gmx_restrict pme)
{
  // unpack params
  const double*  zs = (double*) spline->zs;
  const double*  zx = (double*) spline->theta[0];
  const double*  zy = (double*) spline->theta[1];
  const double*  zz = (double*) spline->theta[2];
  const double* zfx = (double*) spline->dtheta[0];
  const double* zfy = (double*) spline->dtheta[1];
  const double* zfz = (double*) spline->dtheta[2];

  /* ASSUME P=8 const int p = params->P; */
  const int N = params->N;
  const double h=params->h;

  int i,j,idx_zs,m,mm;
  double qm, h3=h*h*h;
  int * idxptr;
  int i0,j0,k0;
  int index_x, index_xy;

  int pny   = pme->pmegrid_ny;
  int pnz   = pme->pmegrid_nz;

  // hold entire zz vector
  __m256d rZZ0, rZZ1;
  __m256d rC, rCX, rCY;
  __m256d rH0, rH1; 
  __m256d rZS0, rZS1;
  __m256d rZFZ0, rZFZ1;
  __m256d rFX, rFY, rFZ;
    
#ifdef CALC_ENERGY
  __m256d rP, rCP;
#endif

  for(m=0; m<N; m++){
    mm  = spline->ind[m];
    qm  = q[mm];
    idxptr = atc->idx[mm];
    i0 = idxptr[XX];
    j0 = idxptr[YY];
    k0 = idxptr[ZZ];

    idx_zs = 0;
    rFX = _mm256_setzero_pd();
    rFY = _mm256_setzero_pd();
    rFZ = _mm256_setzero_pd();
#ifdef CALC_ENERGY
    rP  = _mm256_setzero_pd();
#endif

    /* hoist load of ZZ vector */
    rZZ0 = _mm256_load_pd(zz + m*8     );
    rZZ1 = _mm256_load_pd(zz + m*8 + 4 );

    /* hoist load of ZFZ vector */
    rZFZ0 = _mm256_load_pd(zfz + m*8     );
    rZFZ1 = _mm256_load_pd(zfz + m*8 + 4 );

    if(0){ // H[idx] is 32-aligned so vectorization simple
      for(i = 0; i<8; i++){
	index_x = (i0+i)*pny*pnz;
	real zxi = zx[m*8+i];
	for(j = 0; j<8; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  real zxzy = zxi*zy[m*8+j];
	  rC  = _mm256_set1_pd( zxzy);
	  rCX = _mm256_set1_pd( zxzy * zfx[m*8 + i]);
	  rCY = _mm256_set1_pd( zxzy * zfy[m*8 + j]);
#ifdef CALC_ENERGY
	  rCP  = _mm256_set1_pd( zxzy);
#endif

	  rH0  = _mm256_load_pd( grid + index_xy );
	  rH1  = _mm256_load_pd( grid + index_xy + 4 );
		 
	  rZS0 = _mm256_load_pd( zs + idx_zs     );
	  rZS1 = _mm256_load_pd( zs + idx_zs + 4 );

	  __m256d fac0 = _mm256_mul_pd(rH0, _mm256_mul_pd(rZZ0, rZS0));
	  __m256d fac1 = _mm256_mul_pd(rH1, _mm256_mul_pd(rZZ1, rZS1));
	  __m256d fac = _mm256_add_pd(fac0,fac1);

	  rFX = _mm256_add_pd(rFX,_mm256_mul_pd(fac,rCX));
	  rFY = _mm256_add_pd(rFY,_mm256_mul_pd(fac,rCY));

	  fac0 = _mm256_mul_pd(fac0,rZFZ0);
	  fac1 = _mm256_mul_pd(fac1,rZFZ1);
	  
	  rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(_mm256_add_pd(fac0,fac1),rC));	 

#ifdef CALC_ENERGY
	  rP =_mm256_add_pd(rP,_mm256_mul_pd(fac,rC));
#endif
	  idx_zs +=8;
	}
      }
    }
    else{ // H[idx] not 32-aligned, so use non-aligned loads
      for(i = 0; i<8; i++){
	index_x = (i0+i)*pny*pnz;
	real zxi = zx[m*8+i];
	for(j = 0; j<8; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  real zxzy = zxi*zy[m*8+j];
	  rC  = _mm256_set1_pd( zxzy);
	  rCX = _mm256_set1_pd( zxzy * zfx[m*8 + i]);
	  rCY = _mm256_set1_pd( zxzy * zfy[m*8 + j]);
#ifdef CALC_ENERGY
	  rCP  = _mm256_set1_pd( zxzy);
#endif

	  rH0  = _mm256_loadu_pd( grid + index_xy );
	  rH1  = _mm256_loadu_pd( grid + index_xy + 4 );
		 
	  rZS0 = _mm256_load_pd( zs + idx_zs     );
	  rZS1 = _mm256_load_pd( zs + idx_zs + 4 );

	  __m256d fac0 = _mm256_mul_pd(rH0, _mm256_mul_pd(rZZ0, rZS0));
	  __m256d fac1 = _mm256_mul_pd(rH1, _mm256_mul_pd(rZZ1, rZS1));
	  __m256d fac = _mm256_add_pd(fac0,fac1);

	  rFX = _mm256_add_pd(rFX,_mm256_mul_pd(fac,rCX));
	  rFY = _mm256_add_pd(rFY,_mm256_mul_pd(fac,rCY));

	  fac0 = _mm256_mul_pd(fac0,rZFZ0);
	  fac1 = _mm256_mul_pd(fac1,rZFZ1);
	  
	  rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(_mm256_add_pd(fac0,fac1),rC));	 

#ifdef CALC_ENERGY
	  rP =_mm256_add_pd(rP,_mm256_mul_pd(fac,rC));
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
    force[m][XX] += factor*(reduce(rFX));
    force[m][YY] += factor*(reduce(rFY));
    force[m][ZZ] += factor*(reduce(rFZ));

#ifdef CALC_ENERGY
    st->phi[m] = -scale*h3*reduce(rP);
#endif
  }
}

// -----------------------------------------------------------------------------
static void 
SE_int_split_AVX_P16_gaussian_d(rvec * gmx_restrict force, real * gmx_restrict grid, 
				real * gmx_restrict q,
				splinedata_t          * gmx_restrict spline,
				const SE_params   * gmx_restrict params, 
				real scale, gmx_bool bClearF,
				const pme_atomcomm_t  * gmx_restrict atc,
				const gmx_pme_t       * gmx_restrict pme)
{
  // unpack params
  const double*  zs = (double*) spline->zs;
  const double*  zx = (double*) spline->theta[0];
  const double*  zy = (double*) spline->theta[1];
  const double*  zz = (double*) spline->theta[2];
  const double* zfx = (double*) spline->dtheta[0];
  const double* zfy = (double*) spline->dtheta[1];
  const double* zfz = (double*) spline->dtheta[2];

  /* ASSUME P=16 const int p = params->P; */
  const int N = params->N;
  const double h=params->h;
   
  int i,j,idx,idx_zs,m,mm;
  double qm, h3=h*h*h;
  int * idxptr;
  int i0,j0,k0;
  int index_x, index_xy;

  int pny   = pme->pmegrid_ny;
  int pnz   = pme->pmegrid_nz;

  // hold entire zz vector
  __m256d rZZ0, rZZ1, rZZ2, rZZ3;
  __m256d rC, rCX, rCY;
  __m256d rH0, rH1, rH2, rH3; 
  __m256d rZS0, rZS1, rZS2, rZS3;
  __m256d rZFZ0, rZFZ1, rZFZ2, rZFZ3;
  __m256d rFX, rFY, rFZ;
    
#ifdef CALC_ENERGY
  __m256d rP, rCP;
#endif

  for(m=0; m<N; m++){
    idx = spline->idx[m];
    mm  = spline->ind[m];
    qm = q[mm];
    idxptr = atc->idx[mm];
    i0 = idxptr[XX];
    j0 = idxptr[YY];
    k0 = idxptr[ZZ];
    idx = i0*pny*pnz+j0*pnz+k0;

    _mm_prefetch( (void*) (grid + idx), _MM_HINT_T0 );
    idx_zs = 0;
    _mm_prefetch( (void*) zs, _MM_HINT_T0 );
    rFX = _mm256_setzero_pd();
    rFY = _mm256_setzero_pd();
    rFZ = _mm256_setzero_pd();
#ifdef CALC_ENERGY
    rP  = _mm256_setzero_pd();
#endif


    /* hoist load of ZZ vector */
    rZZ0 = _mm256_load_pd(zz + m*16     );
    rZZ1 = _mm256_load_pd(zz + m*16 + 4 );
    rZZ2 = _mm256_load_pd(zz + m*16 + 8 );
    rZZ3 = _mm256_load_pd(zz + m*16 + 12);

    /* hoist load of ZFZ vector */
    rZFZ0 = _mm256_load_pd(zfz + m*16     );
    rZFZ1 = _mm256_load_pd(zfz + m*16 + 4 );
    rZFZ2 = _mm256_load_pd(zfz + m*16 + 8 );
    rZFZ3 = _mm256_load_pd(zfz + m*16 + 12);

    if(0){ // H[idx] is 32-aligned so vectorization simple
      for(i = 0; i<16; i++){
	index_x = (i0+i)*pny*pnz;
	real zxi = zx[m*16+i];
	for(j = 0; j<16; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  real zxzy = zxi*zy[m*16+j];
	  rC  = _mm256_set1_pd( zxzy );
	  rCX = _mm256_set1_pd( zxzy * zfx[m*16 + i]);
	  rCY = _mm256_set1_pd( zxzy * zfy[m*16 + j]);
#ifdef CALC_ENERGY
	  rCP  = _mm256_set1_pd( zxzy);
#endif

	  rH0  = _mm256_load_pd( grid + index_xy );
	  rH1  = _mm256_load_pd( grid + index_xy + 4 );
	  rH2  = _mm256_load_pd( grid + index_xy + 8 );
	  rH3  = _mm256_load_pd( grid + index_xy + 12);

	  rZS0 = _mm256_load_pd( zs + idx_zs     );
	  rZS1 = _mm256_load_pd( zs + idx_zs + 4 );
	  rZS2 = _mm256_load_pd( zs + idx_zs + 8 );
	  rZS3 = _mm256_load_pd( zs + idx_zs + 12);
		    
	  __m256d fac0 = _mm256_mul_pd(rH0, _mm256_mul_pd(rZZ0, rZS0));
	  __m256d fac1 = _mm256_mul_pd(rH1, _mm256_mul_pd(rZZ1, rZS1));
	  __m256d fac2 = _mm256_mul_pd(rH2, _mm256_mul_pd(rZZ2, rZS2));
	  __m256d fac3 = _mm256_mul_pd(rH3, _mm256_mul_pd(rZZ3, rZS3));

	  __m256d fac = _mm256_add_pd(_mm256_add_pd(fac0,fac1),
				      _mm256_add_pd(fac2,fac3));

	  rFX = _mm256_add_pd(rFX,_mm256_mul_pd(fac,rCX));
	  rFY = _mm256_add_pd(rFY,_mm256_mul_pd(fac,rCY));

	  fac0 = _mm256_mul_pd(fac0,rZFZ0);
	  fac1 = _mm256_mul_pd(fac1,rZFZ1);
	  fac2 = _mm256_mul_pd(fac2,rZFZ2);
	  fac3 = _mm256_mul_pd(fac3,rZFZ3);

	  __m256d fact =_mm256_add_pd( _mm256_add_pd(fac0,fac1), 
				       _mm256_add_pd(fac2,fac3) );

	  rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(fact,rC));

#ifdef CALC_ENERGY
	  rP =_mm256_add_pd(rP,_mm256_mul_pd(fac,rC));
#endif

	  idx_zs +=16;
	}
      }
    }
    else{ // H[idx] not 32-aligned, so use non-aligned loads
      for(i = 0; i<16; i++){
	index_x = (i0+i)*pny*pnz;
	real zxi = zx[m*16+i];
	for(j = 0; j<16; j++){
	  index_xy = index_x + (j0+j)*pnz + k0;
	  real zxzy = zxi*zy[m*16+j];
	  rC  = _mm256_set1_pd( zxzy );
	  rCX = _mm256_set1_pd( zxzy * zfx[m*16 + i]);
	  rCY = _mm256_set1_pd( zxzy * zfy[m*16 + j]);
#ifdef CALC_ENERGY
	  rCP  = _mm256_set1_pd( zxzy);
#endif

	  rH0  = _mm256_loadu_pd( grid + index_xy );
	  rH1  = _mm256_loadu_pd( grid + index_xy + 4 );
	  rH2  = _mm256_loadu_pd( grid + index_xy + 8 );
	  rH3  = _mm256_loadu_pd( grid + index_xy + 12);

	  rZS0 = _mm256_load_pd( zs + idx_zs     );
	  rZS1 = _mm256_load_pd( zs + idx_zs + 4 );
	  rZS2 = _mm256_load_pd( zs + idx_zs + 8 );
	  rZS3 = _mm256_load_pd( zs + idx_zs + 12);
		    
	  __m256d fac0 = _mm256_mul_pd(rH0, _mm256_mul_pd(rZZ0, rZS0));
	  __m256d fac1 = _mm256_mul_pd(rH1, _mm256_mul_pd(rZZ1, rZS1));
	  __m256d fac2 = _mm256_mul_pd(rH2, _mm256_mul_pd(rZZ2, rZS2));
	  __m256d fac3 = _mm256_mul_pd(rH3, _mm256_mul_pd(rZZ3, rZS3));

	  __m256d fac = _mm256_add_pd(_mm256_add_pd(fac0,fac1),
				      _mm256_add_pd(fac2,fac3));

	  rFX = _mm256_add_pd(rFX,_mm256_mul_pd(fac,rCX));
	  rFY = _mm256_add_pd(rFY,_mm256_mul_pd(fac,rCY));

	  fac0 = _mm256_mul_pd(fac0,rZFZ0);
	  fac1 = _mm256_mul_pd(fac1,rZFZ1);
	  fac2 = _mm256_mul_pd(fac2,rZFZ2);
	  fac3 = _mm256_mul_pd(fac3,rZFZ3);

	  __m256d fact =_mm256_add_pd( _mm256_add_pd(fac0,fac1), 
				       _mm256_add_pd(fac2,fac3) );

	  rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(fact,rC));

#ifdef CALC_ENERGY
	  rP =_mm256_add_pd(rP,_mm256_mul_pd(fac,rC));
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
    force[m][XX] += factor*(reduce(rFX));
    force[m][YY] += factor*(reduce(rFY));
    force[m][ZZ] += factor*(reduce(rFZ));

#ifdef CALC_ENERGY
    st->phi[m] = -scale*h3*reduce(rP);
#endif
  }

}

/* =============================================================================
 * =============================== KAISER ROUTINES =============================
 * =============================================================================
 */
static
void SE_int_split_AVX_P4_kaiser_d(rvec * gmx_restrict force, real * gmx_restrict grid, 
				  real * gmx_restrict q,
				  splinedata_t          * gmx_restrict spline,
				  const SE_params   * gmx_restrict params, 
				  real scale, gmx_bool bClearF,
				  const pme_atomcomm_t  * gmx_restrict atc,
				  const gmx_pme_t       * gmx_restrict pme)
{  
  // unpack params
  const double* zx = (double*) spline->theta[0];
  const double* zy = (double*) spline->theta[1];
  const double* zz = (double*) spline->theta[2];
  const double* zfx =(double*) spline->dtheta[0];
  const double* zfy =(double*) spline->dtheta[1];
  const double* zfz =(double*) spline->dtheta[2];

  const int    p = params->P;
  const int    N = params->N;
  const double h = params->h;

  int i,j, m;
  double qm, h3=h*h*h;
  int * idxptr;
  int i0,j0,k0;
  int index_x, index_xy;

  int pny   = pme->pmegrid_ny;
  int pnz   = pme->pmegrid_nz;

  __m256d rH0, rZZ0, rZFZ0;
  __m256d rC, rCX, rCY;
  __m256d rFX, rFY, rFZ;
  
    for(int mm=0; mm<N; mm++)
      {
	m  = spline->ind[mm];
	qm  = q[m];
	idxptr = atc->idx[m];
	i0 = idxptr[XX];
	j0 = idxptr[YY];
	k0 = idxptr[ZZ];

	rFX = _mm256_setzero_pd();
	rFY = _mm256_setzero_pd();
	rFZ = _mm256_setzero_pd();
	
	int m4 = m*4;
        // hoist load of ZZ and ZFZ vector //
        rZZ0 = _mm256_load_pd( zz + m4  );
	rZFZ0= _mm256_load_pd( zfz+ m*p );
	
	for(i = 0; i<4; i++)
	  {
	    index_x = (i0+i)*pny*pnz;
	    real zxi = zx[m4+i];
	    for(j = 0; j<4; j++)
	      {
		index_xy = index_x + (j0+j)*pnz + k0;
		real zxzy = zxi*zy[m4+j];
		rC  = _mm256_set1_pd( zxzy );
		rCX = _mm256_set1_pd( zxzy*zfx[m4+i] );
		rCY = _mm256_set1_pd( zxzy*zfy[m4+j] );
		
		rH0  = _mm256_loadu_pd( grid + index_xy );

		__m256d fac = _mm256_mul_pd(rH0, rZZ0);

#ifdef AVX_FMA
		rFX = _mm256_fmadd_pd(fac,rCX,rFX);
		rFY = _mm256_fmadd_pd(fac,rCY,rFY);
		rFZ = _mm256_fmadd_pd(fac,_mm256_mul_pd(rC,rZFZ0),rFZ);
#else
		rFX = _mm256_add_pd(_mm256_mul_pd(fac,rCX), rFX);
		rFY = _mm256_add_pd(_mm256_mul_pd(fac,rCY), rFY);
		rFZ = _mm256_add_pd(_mm256_mul_pd(fac,_mm256_mul_pd(rC,rZFZ0)), rFZ);
#endif
	      }
        }
	if(bClearF){
	  force[m][XX] = 0;
	  force[m][YY] = 0;
	  force[m][ZZ] = 0;
	}
	
	double factor = -qm*scale*h3;
	force[m][XX] += factor*(reduce(rFX));
	force[m][YY] += factor*(reduce(rFY));
	force[m][ZZ] += factor*(reduce(rFZ));
      }
}

static
void SE_int_split_AVX_P6_kaiser_d(rvec * gmx_restrict force, real * gmx_restrict grid, 
				  real * gmx_restrict q,
				  splinedata_t          * gmx_restrict spline,
				  const SE_params   * gmx_restrict params, 
				  real scale, gmx_bool bClearF,
				  const pme_atomcomm_t  * gmx_restrict atc,
				  const gmx_pme_t       * gmx_restrict pme)
{  
    // unpack params
  const double*  zx = (double*) spline->theta[0];
  const double*  zy = (double*) spline->theta[1];
  const double*  zz = (double*) spline->theta[2];
  const double* zfx = (double*) spline->dtheta[0];
  const double* zfy = (double*) spline->dtheta[1];
  const double* zfz = (double*) spline->dtheta[2];

  /* ASSUME P=8 const int p = params->P; */
  const int N = params->N;
  const double h=params->h;

  int i,j, m;
  double qm, h3=h*h*h;
  int * idxptr;
  int i0,j0,k0;
  int index_x, index_xy;

  int pny   = pme->pmegrid_ny;
  int pnz   = pme->pmegrid_nz;

  // hold entire zz vector
  __m256d rZZ0;
  __m256d rC0, rCX0, rCY0;
  __m256d rH0; 
  __m256d rZFZ0;
  __m256d rFX0, rFY0, rFZ0;

  __m128d rZZ1; 
  __m128d rC1, rCX1, rCY1;
  __m128d rH1; 
  __m128d rZFZ1;
  __m128d rFX1, rFY1, rFZ1;
  
  for(int mm=0; mm<N; mm++)
    {
      m  = spline->ind[mm];
      qm  = q[m];
      idxptr = atc->idx[m];
      i0 = idxptr[XX];
      j0 = idxptr[YY];
      k0 = idxptr[ZZ];
      
      rFX0 = _mm256_setzero_pd();
      rFY0 = _mm256_setzero_pd();
      rFZ0 = _mm256_setzero_pd();

      rFX1 = _mm_setzero_pd();
      rFY1 = _mm_setzero_pd();
      rFZ1 = _mm_setzero_pd();

      int m6 = m*6;
      
      // hoist load of ZZ vector //
      rZZ0 = _mm256_loadu_pd(zz + m6     );
      rZZ1 = _mm_loadu_pd(   zz + m6 + 4 );

      /* hoist load of ZFZ vector */
      rZFZ0 = _mm256_loadu_pd(zfz + m*6     );
      rZFZ1 = _mm_loadu_pd(   zfz + m*6 + 4 );
      
      for(i = 0; i<6; i++)
	{
	  index_x = (i0+i)*pny*pnz;
	  real zxi = zx[m6+i];
	  for(j = 0; j<6; j++)
	    {
	      index_xy = index_x + (j0+j)*pnz + k0;
	      real zxzy = zxi*zy[m6+j];
	      rC0  = _mm256_set1_pd( zxzy);
	      rCX0 = _mm256_set1_pd( zxzy * zfx[m6 + i]);
	      rCY0 = _mm256_set1_pd( zxzy * zfy[m6 + j]);

	      rC1  = _mm_set1_pd( zxzy);
	      rCX1 = _mm_set1_pd( zxzy * zfx[m6 + i]);
	      rCY1 = _mm_set1_pd( zxzy * zfy[m6 + j]);
	      
	      rH0  = _mm256_loadu_pd( grid + index_xy    );
	      rH1  = _mm_loadu_pd(    grid + index_xy + 4);
	      
	      __m256d fac0 = _mm256_mul_pd(rH0, rZZ0);
	      __m128d fac1 = _mm_mul_pd(rH1, rZZ1);

#ifdef AVX_FMA
	      rFX0 = _mm256_fmadd_pd(fac0,rCX0,rFX0);
	      rFY0 = _mm256_fmadd_pd(fac0,rCY0,rFY0);

	      rFX1 = _mm_fmadd_pd(fac1,rCX1,rFX1);
	      rFY1 = _mm_fmadd_pd(fac1,rCY1,rFY1);
#else
	      rFX0 = _mm256_add_pd(_mm256_mul_pd(fac0,rCX0),rFX0);
	      rFY0 = _mm256_add_pd(_mm256_mul_pd(fac0,rCY0),rFY0);

	      rFX1 = _mm_add_pd(_mm_mul_pd(fac1,rCX1),rFX1);
	      rFY1 = _mm_add_pd(_mm_mul_pd(fac1,rCY1),rFY1);
#endif
	      fac0 = _mm256_mul_pd(fac0,rZFZ0);
	      fac1 = _mm_mul_pd(fac1,rZFZ1);
#ifdef AVX_FMA
	      rFZ0 = _mm256_fmadd_pd(fac0,rC0,rFZ0);
	      rFZ1 = _mm_fmadd_pd(fac1,rC1,rFZ1);
#else
	      rFZ0 = _mm256_add_pd(_mm256_mul_pd(fac0,rC0),rFZ0);
	      rFZ1 = _mm_add_pd(_mm_mul_pd(fac1,rC1),rFZ1);
#endif
	    }
        }
      if(bClearF){
	force[m][XX] = 0;
	force[m][YY] = 0;
	force[m][ZZ] = 0;
      }

      double factor = -qm*scale*h3;
      force[m][XX] += factor*(reduce(rFX0)+reduce_sse(rFX1));
      force[m][YY] += factor*(reduce(rFY0)+reduce_sse(rFY1));
      force[m][ZZ] += factor*(reduce(rFZ0)+reduce_sse(rFZ1));
    }
}


static
void SE_int_split_AVX_P8_kaiser_d(rvec * gmx_restrict force, real * gmx_restrict grid, 
				  real * gmx_restrict q,
				  splinedata_t          * gmx_restrict spline,
				  const SE_params   * gmx_restrict params, 
				  real scale, gmx_bool bClearF,
				  const pme_atomcomm_t  * gmx_restrict atc,
				  const gmx_pme_t       * gmx_restrict pme)
{  
    // unpack params
  const double*  zx = (double*) spline->theta[0];
  const double*  zy = (double*) spline->theta[1];
  const double*  zz = (double*) spline->theta[2];
  const double* zfx = (double*) spline->dtheta[0];
  const double* zfy = (double*) spline->dtheta[1];
  const double* zfz = (double*) spline->dtheta[2];

  /* ASSUME P=8 const int p = params->P; */
  const int N = params->N;
  const double h=params->h;

  int i,j, m;
  double qm, h3=h*h*h;
  int * idxptr;
  int i0,j0,k0;
  int index_x, index_xy;

  int pny   = pme->pmegrid_ny;
  int pnz   = pme->pmegrid_nz;

  // hold entire zz vector
  __m256d rZZ0, rZZ1;
  __m256d rC, rCX, rCY;
  __m256d rH0, rH1; 
  __m256d rZFZ0, rZFZ1;
  __m256d rFX, rFY, rFZ;
  
    for(int mm=0; mm<N; mm++)
    {
      m  = spline->ind[mm];
      qm  = q[m];
      idxptr = atc->idx[m];
      i0 = idxptr[XX];
      j0 = idxptr[YY];
      k0 = idxptr[ZZ];

      rFX = _mm256_setzero_pd();
      rFY = _mm256_setzero_pd();
      rFZ = _mm256_setzero_pd();

      int m8 = m*8;
      
      // hoist load of ZZ vector //
      rZZ0 = _mm256_load_pd(zz + m8     );
      rZZ1 = _mm256_load_pd(zz + m8 + 4 );

      /* hoist load of ZFZ vector */
      rZFZ0 = _mm256_load_pd(zfz + m8     );
      rZFZ1 = _mm256_load_pd(zfz + m8 + 4 );
      
      for(i = 0; i<8; i++)
	{
	  index_x = (i0+i)*pny*pnz;
	  real zxi = zx[m8+i];
	  for(j = 0; j<8; j++)
	    {
	      index_xy = index_x + (j0+j)*pnz + k0;
	      real zxzy = zxi*zy[m8+j];
	      rC  = _mm256_set1_pd( zxzy);
	      rCX = _mm256_set1_pd( zxzy * zfx[m8 + i]);
	      rCY = _mm256_set1_pd( zxzy * zfy[m8 + j]);
	      
	      rH0  = _mm256_loadu_pd( grid + index_xy );
	      rH1  = _mm256_loadu_pd( grid + index_xy + 4 );

	      __m256d fac0 = _mm256_mul_pd(rH0, rZZ0);
	      __m256d fac1 = _mm256_mul_pd(rH1, rZZ1);
	      __m256d fac = _mm256_add_pd(fac0,fac1);

#ifdef AVX_FMA
	      rFX = _mm256_fmadd_pd(fac,rCX,rFX);
	      rFY = _mm256_fmadd_pd(fac,rCY,rFY);	      
#else
	      rFX = _mm256_add_pd(_mm256_mul_pd(fac,rCX),rFX);
	      rFY = _mm256_add_pd(_mm256_mul_pd(fac,rCY),rFY);
#endif
	      fac0 = _mm256_mul_pd(fac0,rZFZ0);
	      fac1 = _mm256_mul_pd(fac1,rZFZ1);
#ifdef AVX_FMA	      
	      rFZ = _mm256_fmadd_pd(_mm256_add_pd(fac0,fac1),rC,rFZ);
#else
	      rFZ = _mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(fac0,fac1),rC),rFZ);
#endif
	    }
        }
      if(bClearF){
	force[m][XX] = 0;
	force[m][YY] = 0;
	force[m][ZZ] = 0;
      }

      double factor = -qm*scale*h3;
      force[m][XX] += factor*(reduce(rFX));
      force[m][YY] += factor*(reduce(rFY));
      force[m][ZZ] += factor*(reduce(rFZ));
    }
}

static
void SE_int_split_AVX_P16_kaiser_d(rvec * gmx_restrict force, real * gmx_restrict grid, 
				   real * gmx_restrict q,
				   splinedata_t          * gmx_restrict spline,
				   const SE_params   * gmx_restrict params, 
				   real scale, gmx_bool bClearF,
				   const pme_atomcomm_t  * gmx_restrict atc,
				   const gmx_pme_t       * gmx_restrict pme)
{
    // unpack params
  const double*  zx = (double*) spline->theta[0];
  const double*  zy = (double*) spline->theta[1];
  const double*  zz = (double*) spline->theta[2];
  const double* zfx = (double*) spline->dtheta[0];
  const double* zfy = (double*) spline->dtheta[1];
  const double* zfz = (double*) spline->dtheta[2];

  /* ASSUME P=16 const int p = params->P; */
  const int N = params->N;
  const double h=params->h;
   
  int i,j,m;
  double qm, h3=h*h*h;
  int * idxptr;
  int i0,j0,k0;
  int index_x, index_xy;

  int pny   = pme->pmegrid_ny;
  int pnz   = pme->pmegrid_nz;

  // hold entire zz vector
  __m256d rZZ0, rZZ1, rZZ2, rZZ3;
  __m256d rC, rCX, rCY;
  __m256d rH0, rH1, rH2, rH3; 
  __m256d rZFZ0, rZFZ1, rZFZ2, rZFZ3;
  __m256d rFX, rFY, rFZ;
  
  for(int mm=0; mm<N; mm++)
    {
      m  = spline->ind[mm];
      qm = q[m];
      idxptr = atc->idx[m];
      i0 = idxptr[XX];
      j0 = idxptr[YY];
      k0 = idxptr[ZZ];

      int m16 = m*16;

      rFX = _mm256_setzero_pd();
      rFY = _mm256_setzero_pd();
      rFZ = _mm256_setzero_pd();
      
      /* hoist load of ZZ vector */
      rZZ0 = _mm256_load_pd(zz + m16     );
      rZZ1 = _mm256_load_pd(zz + m16 + 4 );
      rZZ2 = _mm256_load_pd(zz + m16 + 8 );
      rZZ3 = _mm256_load_pd(zz + m16 + 12);

      /* hoist load of ZFZ vector */
      rZFZ0 = _mm256_load_pd(zfz + m16     );
      rZFZ1 = _mm256_load_pd(zfz + m16 + 4 );
      rZFZ2 = _mm256_load_pd(zfz + m16 + 8 );
      rZFZ3 = _mm256_load_pd(zfz + m16 + 12);
      
      for(i = 0; i<16; i++)
	{
	  index_x = (i0+i)*pny*pnz;
	  real zxi = zx[m16+i];
	  for(j = 0; j<16; j++)
	    {
	      index_xy = index_x + (j0+j)*pnz + k0;
	      real zxzy = zxi*zy[m16+j];
	      rC  = _mm256_set1_pd( zxzy );
	      rCX = _mm256_set1_pd( zxzy * zfx[m16 + i]);
	      rCY = _mm256_set1_pd( zxzy * zfy[m16 + j]);

	      rH0  = _mm256_loadu_pd( grid + index_xy );
	      rH1  = _mm256_loadu_pd( grid + index_xy + 4 );
	      rH2  = _mm256_loadu_pd( grid + index_xy + 8 );
	      rH3  = _mm256_loadu_pd( grid + index_xy + 12);

	      __m256d fac0 = _mm256_mul_pd(rH0,rZZ0);
	      __m256d fac1 = _mm256_mul_pd(rH1,rZZ1);
	      __m256d fac2 = _mm256_mul_pd(rH2,rZZ2);
	      __m256d fac3 = _mm256_mul_pd(rH3,rZZ3);

	      __m256d fac = _mm256_add_pd(_mm256_add_pd(fac0,fac1),
					  _mm256_add_pd(fac2,fac3));

#ifdef AVX_FMA
	      rFX = _mm256_fmadd_pd(fac,rCX,rFX);
	      rFY = _mm256_fmadd_pd(fac,rCY,rFY);
#else
	      rFX = _mm256_add_pd(_mm256_mul_pd(fac,rCX),rFX);
	      rFY = _mm256_add_pd(_mm256_mul_pd(fac,rCY),rFY);	      
#endif
	      
	      fac0 = _mm256_mul_pd(fac0,rZFZ0);
	      fac1 = _mm256_mul_pd(fac1,rZFZ1);
	      fac2 = _mm256_mul_pd(fac2,rZFZ2);
	      fac3 = _mm256_mul_pd(fac3,rZFZ3);
	      
	      __m256d fact =_mm256_add_pd( _mm256_add_pd(fac0,fac1), 
					   _mm256_add_pd(fac2,fac3) );
#ifdef AVX_FMA
	      rFZ = _mm256_fmadd_pd(fact,rC,rFZ);
#else
	      rFZ = _mm256_add_pd(_mm256_mul_pd(fact,rC),rFZ);	      
#endif
	    }
	}
      if(bClearF){
	force[m][XX] = 0;
	force[m][YY] = 0;
	force[m][ZZ] = 0;
      }

      double factor = -qm*scale*h3;
      force[m][XX] += factor*(reduce(rFX));
      force[m][YY] += factor*(reduce(rFY));
      force[m][ZZ] += factor*(reduce(rFZ));
    }
}
// -----------------------------------------------------------------------------
static void 
SE_int_split_AVX_dispatch_d(rvec* force, real* grid, real* q, 
			    splinedata_t *spline,
			    const SE_params* params, real scale,
			    gmx_bool bClearF,
			    const pme_atomcomm_t *atc,
			    const gmx_pme_t *pme,
			    const int se_set)
{
  const int p = params->P;
  const int incrj = params->dims[2]; // middle increment
  const int incri = params->npdims[2]*(params->dims[1]);// outer increment
  
  if(se_set==1)
    {
      // if P, incri or increments are not divisible by 4, fall back on vanilla
      if(isnot_div_by_4(p) || isnot_div_by_4(incri) || isnot_div_by_4(incrj) ){
	__DISPATCHER_MSG("[FGG INT AVX DOUBLE] AVX Abort (PARAMS)\n");
	SE_int_split_SSE_dispatch_d(force, grid, q, spline, params,
				    scale, bClearF, atc, pme, se_set);
	return;
      }
    }
  else
    {
      if(isnot_div_by_4(incri) || isnot_div_by_4(incrj) ){
	__DISPATCHER_MSG("[FKG INT AVX DOUBLE] AVX Abort (PARAMS)\n");
	SE_int_split_SSE_dispatch_d(force, grid, q, spline, params,
				    scale, bClearF, atc, pme, se_set);
	return;
      }
    }

  if(se_set==1) // Gaussian 
    {
      // otherwise the preconditions for AVX codes are satisfied. 
      if(p==8){
	// specific for p=8
	__DISPATCHER_MSG("[FGG INT AVX DOUBLE] P=8\n");
	SE_int_split_AVX_P8_gaussian_d(force, grid, q, spline, params,
				       scale, bClearF, atc, pme);
      }
      else if(p==16){
	// specific for p=16
	__DISPATCHER_MSG("[FGG INT AVX DOUBLE] P=16\n");
	SE_int_split_AVX_P16_gaussian_d(force, grid, q, spline, params,
					scale, bClearF, atc, pme);
      }
      else if(p%8==0){
	// for p divisible by 8
	__DISPATCHER_MSG("[FGG INT AVX DOUBLE] P unroll 8\n");
	SE_int_split_AVX_u8_gaussian_d(force, grid, q, spline, params,
				       scale, bClearF, atc, pme); 
      }
      else if(p%4==0){
	// vanilla AVX code (p divisible by 4)
	__DISPATCHER_MSG("[FGG INT AVX DOUBLE] P unroll 4\n");
	SE_int_split_AVX_gaussian_d(force, grid, q, spline, params,
				       scale, bClearF, atc, pme);
      }
      else {
	// vanilla SSE code (any even p)
	__DISPATCHER_MSG("[FGG INT AVX DOUBLE] AVX Abort\n");
	SE_int_split_SSE_dispatch_d(force, grid, q, spline, params,
				    scale, bClearF, atc, pme, se_set);
      }
    } //se_set=1 Gaussian
  else // se_set=2 Kaiser
    {
      // otherwise the preconditions for AVX codes are satisfied.
      if(p==4){
      	// specific for p=4
      	__DISPATCHER_MSG("[FKG INT AVX DOUBLE] P unroll 4\n");
      	SE_int_split_AVX_P4_kaiser_d(force, grid, q, spline, params,
				     scale, bClearF, atc, pme);
      }
      else if(p==6) {
      	// specific for p=6
      	__DISPATCHER_MSG("[FKG INT AVX DOUBLE] P=6\n");
      	SE_int_split_AVX_P6_kaiser_d(force, grid, q, spline, params,
				     scale, bClearF, atc, pme);
      }
      else if(p==8) {
      	// specific for p=8
      	__DISPATCHER_MSG("[FKG INT AVX DOUBLE] P=8\n");
      	SE_int_split_AVX_P8_kaiser_d(force, grid, q, spline, params,
				     scale, bClearF, atc, pme);
      }
      else if(p==16){
      	// specific for p=16
      	__DISPATCHER_MSG("[FKG INT AVX DOUBLE] P=16\n");
      	SE_int_split_AVX_P16_kaiser_d(force, grid, q, spline, params,
      				      scale, bClearF, atc, pme);
      }
      else {
      	// vanilla SSE code (any even p)
      	__DISPATCHER_MSG("[FKG INT AVX DOUBLE] AVX Abort\n");
      	SE_int_split_SSE_dispatch_d(force, grid, q, spline, params,
				    scale, bClearF, atc, pme, se_set);
      }
    }
}

#endif //GMX_SIMD_X86_AVX_256
#endif //GMX_DOUBLE
#endif //_SE_INT_AVX_256_DOUBLE_
