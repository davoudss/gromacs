#ifndef _SE_INT_AVX_256_DOUBLE_
#define _SE_INT_AVX_256_DOUBLE_


/* SE AVX 256 double integration */
#if GMX_DOUBLE==1
#if GMX_SIMD_X86_AVX_256
#include "se.h"


static void
SE_int_split_AVX_d(rvec *force, real *grid, real *q,
		   splinedata_t *spline,
		   const SE_FGG_params *params, real scale,
		   gmx_bool bClearF)
{
  // unpack params
  const double*  H = (double*) grid;
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

  int i,j,k,m,idx,idx_zs,idx_zz, mm;
  double qm, h3=h*h*h;
  double sx[4] MEM_ALIGNED;
  double sy[4] MEM_ALIGNED;
  double sz[4] MEM_ALIGNED;

  __m256d rH0, rZZ0, rZS0,rZFZ0;
  __m256d rC, rCX, rCY;
  __m256d rFX, rFY, rFZ;
#ifdef CALC_ENERGY
  double s[4]  MEM_ALIGNED;
  __m256d rP, rCP;
#endif

  const int incrj = params->npdims[2]-p;
  const int incri = params->npdims[2]*(params->npdims[1]-p);

  for(m=0; m<N; m++){
    idx = spline->idx[m];
    mm  = spline->ind[m];
    qm  = q[mm];
    
    idx_zs = 0;
    rFX = _mm256_setzero_pd();
    rFY = _mm256_setzero_pd();
    rFZ = _mm256_setzero_pd();
#ifdef CALC_ENERGY
    rP = _mm256_setzero_pd();
#endif
    
    if(idx%4==0){ // H[idx] is 32-aligned so vectorization simple
      for(i = 0; i<p; i++){
	for(j = 0; j<p; j++){
	  rC  = _mm256_set1_pd( zx[m*p+i]*zy[m*p+j] );
	  rCX = _mm256_set1_pd(zx[m*p+i]*zy[m*p+j]*zfx[m*p+i]);
	  rCY = _mm256_set1_pd(zx[m*p+i]*zy[m*p+j]*zfy[m*p+j]);
#ifdef CALC_ENERGY
	  rCP  = _mm256_set1_pd( zx[m*p+i]*zy[m*p+j]);
#endif

	  idx_zz=m*p;
	  for(k = 0; k<p; k+=4){
	    rZFZ0= _mm256_load_pd( zfz+ m*p+k );
	    rH0  = _mm256_load_pd( H  + idx );
	    rZZ0 = _mm256_load_pd( zz + idx_zz);
	    rZS0 = _mm256_load_pd( zs + idx_zs);
	    rFX = _mm256_add_pd(rFX,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCX),rZS0)));
	    rFY = _mm256_add_pd(rFY,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCY),rZS0)));
	    rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)),rZFZ0));

#ifdef CALC_ENERGY
	    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCP),rZS0)));
#endif			
	    idx+=4; 
	    idx_zs+=4; 
	    idx_zz+=4;
	  }
	  idx += incrj;
	}
	idx += incri;
      }
    }
    else{ // H[idx] not 16-aligned, so use non-aligned loads
      for(i = 0; i<p; i++){
	for(j = 0; j<p; j++){
	  rC  = _mm256_set1_pd( zx[m*p+i]*zy[m*p+j] );
	  rCX = _mm256_set1_pd(zx[m*p+i]*zy[m*p+j]*zfx[m*p+i]);
	  rCY = _mm256_set1_pd(zx[m*p+i]*zy[m*p+j]*zfy[m*p+j]);
#ifdef CALC_ENERGY
	  rCP = _mm256_set1_pd( zx[m*p+i]*zy[m*p+j] );
#endif
	  idx_zz=m*p;
	  for(k = 0; k<p; k+=4){
	    rZFZ0= _mm256_load_pd( zfz + m*p+k );
	    rH0  = _mm256_loadu_pd( H+idx );
	    rZZ0 = _mm256_load_pd( zz + idx_zz);
	    rZS0 = _mm256_load_pd( zs + idx_zs);
	    rFX = _mm256_add_pd(rFX,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCX),rZS0)));
	    rFY = _mm256_add_pd(rFY,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCY),rZS0)));
	    rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)),rZFZ0));

#ifdef CALC_ENERGY
	    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCP),rZS0)));
#endif			
	    idx+=4; 
	    idx_zs+=4; 
	    idx_zz+=4;
	  }
	  idx += incrj;
	}
	idx += incri;
      }

    }
    _mm256_store_pd(sx,rFX);
    _mm256_store_pd(sy,rFY);
    _mm256_store_pd(sz,rFZ);

    if(bClearF){
      force[m][XX] = 0;
      force[m][YY] = 0;
      force[m][ZZ] = 0;
    }

    force[m][XX] += -qm*scale*h3*(sx[0]+sx[1]+sx[2]+sx[3]);
    force[m][YY] += -qm*scale*h3*(sy[0]+sy[1]+sy[2]+sy[3]);
    force[m][ZZ] += -qm*scale*h3*(sz[0]+sz[1]+sz[2]+sz[3]);

#ifdef CALC_ENERGY
    _mm256_store_pd(s,rP);
    st->phi[m] = -scale*h3*(s[0]+s[1]+s[2]+s[3]);
#endif
  }
}

// -----------------------------------------------------------------------------
static void 
SE_int_split_AVX_u8_d(rvec* force, real* grid, real* q,  
		      splinedata_t *spline,
		      const SE_FGG_params* params, real scale,
		      gmx_bool bClearF)
{
  // unpack params
  const double*   H = (double*) grid;
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
  double sx[4] MEM_ALIGNED;
  double sy[4] MEM_ALIGNED;
  double sz[4] MEM_ALIGNED;

  __m256d rH0, rZZ0, rZS0, rZFZ0;
  __m256d rH1, rZZ1, rZS1, rZFZ1;
  __m256d rFX, rFY, rFZ;
  __m256d  rC, rCX, rCY;

#ifdef CALC_ENERGY
  double s[4]  MEM_ALIGNED;
  __m256d rP, rCP;
#endif

  const int incrj = params->npdims[2]-p;
  const int incri = params->npdims[2]*(params->npdims[1]-p);

  for(m=0; m<N; m++){
    idx = spline->idx[m];
    mm = spline->ind[m];
    qm = q[mm];

    _mm_prefetch( (void*) (H+idx), _MM_HINT_T0);
    idx_zs = 0;
    _mm_prefetch( (void*) zs, _MM_HINT_T0);
    rFX = _mm256_setzero_pd();
    rFY = _mm256_setzero_pd();
    rFZ = _mm256_setzero_pd();
#ifdef CALC_ENERGY
    rP  = _mm256_setzero_pd();
#endif

    if(idx%4==0){ // H[idx] is 32-aligned so vectorization simple
      for(i = 0; i<p; i++){
	for(j = 0; j<p; j++){
	  rC  = _mm256_set1_pd( zx[m*p+i]*zy[m*p+j]);
	  rCX = _mm256_set1_pd(zx[m*p+i]*zy[m*p+j]*zfx[m*p+i]);
	  rCY = _mm256_set1_pd(zx[m*p+i]*zy[m*p+j]*zfy[m*p+j]);
#ifdef CALC_ENERGY
	  rCP = _mm256_set1_pd( zx[m*p+i]*zy[m*p+j] );
#endif

	  idx_zz=m*p;

	  for(k = 0; k<p; k+=8){
	    rH0  = _mm256_load_pd( H+idx    );
	    rH1  = _mm256_load_pd( H+idx + 4);

	    rZZ0 = _mm256_load_pd( zz + idx_zz    );
	    rZZ1 = _mm256_load_pd( zz + idx_zz + 4);

	    rZS0 = _mm256_load_pd( zs + idx_zs    );
	    rZS1 = _mm256_load_pd( zs + idx_zs + 4);

	    rZFZ0 = _mm256_load_pd(zfz+ idx_zz    );
	    rZFZ1 = _mm256_load_pd(zfz+ idx_zz + 4);

	    rFX = _mm256_add_pd(rFX,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCX),rZS0)));
	    rFX = _mm256_add_pd(rFX,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCX),rZS1)));

	    rFY = _mm256_add_pd(rFY,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCY),rZS0)));
	    rFY = _mm256_add_pd(rFY,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCY),rZS1)));
			
	    rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
	    rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ1,rZZ1),rC),rZS1)));

#ifdef CALC_ENERGY
	    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCP),rZS0)));
	    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCP),rZS1)));
#endif


	    idx+=8; 
	    idx_zs+=8; 
	    idx_zz+=8;
	  }
	  idx += incrj;
	}
	idx += incri;
      }
    }
    else{ // H[idx] not 32-aligned, so use non-aligned load from H
      for(i = 0; i<p; i++){
	for(j = 0; j<p; j++){
	  rC  = _mm256_set1_pd( zx[m*p+i]*zy[m*p+j] );
	  rCX = _mm256_set1_pd(zx[m*p+i]*zy[m*p+j]*zfx[m*p+i]);
	  rCY = _mm256_set1_pd(zx[m*p+i]*zy[m*p+j]*zfy[m*p+j]);
#ifdef CALC_ENERGY
	  rCP = _mm256_set1_pd( zx[m*p+i]*zy[m*p+j] );
#endif

	  idx_zz=m*p;
	  for(k = 0; k<p; k+=8){
	    rH0  = _mm256_loadu_pd( H+idx    );
	    rH1  = _mm256_loadu_pd( H+idx + 4);
		
	    rZZ0 = _mm256_load_pd( zz + idx_zz    );
	    rZZ1 = _mm256_load_pd( zz + idx_zz + 4);

	    rZS0 = _mm256_load_pd( zs + idx_zs    );
	    rZS1 = _mm256_load_pd( zs + idx_zs + 4);

	    rZFZ0 = _mm256_load_pd(zfz+ idx_zz    );
	    rZFZ1 = _mm256_load_pd(zfz+ idx_zz + 4);

	    rFX = _mm256_add_pd(rFX,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCX),rZS0)));
	    rFX = _mm256_add_pd(rFX,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCX),rZS1)));

	    rFY = _mm256_add_pd(rFY,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCY),rZS0)));
	    rFY = _mm256_add_pd(rFY,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCY),rZS1)));
			
	    rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
	    rFZ = _mm256_add_pd(rFZ,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ1,rZZ1),rC),rZS1)));

#ifdef CALC_ENERGY
	    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCP),rZS0)));
	    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCP),rZS1)));
#endif

	    idx+=8; 
	    idx_zs+=8; 
	    idx_zz+=8;
	  }
	  idx += incrj;
	}
	idx += incri;
      }
    }

    // done accumulating
    _mm256_store_pd(sx,rFX);
    _mm256_store_pd(sy,rFY);
    _mm256_store_pd(sz,rFZ);

    if(bClearF){
      force[m][XX] = 0;
      force[m][YY] = 0;
      force[m][ZZ] = 0;
    }

    force[m][XX] += -qm*scale*h3*(sx[0]+sx[1]+sx[2]+sx[3]);
    force[m][YY] += -qm*scale*h3*(sy[0]+sy[1]+sy[2]+sy[3]);
    force[m][ZZ] += -qm*scale*h3*(sz[0]+sz[1]+sz[2]+sz[3]);

#ifdef CALC_ENERGY
    _mm256_store_pd(s,rP);
    st->phi[m] = -scale*h3*(s[0]+s[1]+s[2]+s[3]);
#endif
  }
}

// -----------------------------------------------------------------------------
static void
SE_int_split_AVX_P8_d(rvec* force, real* grid, real* q,  
		      splinedata_t *spline,
		      const SE_FGG_params* params, real scale,
		      gmx_bool bClearF)
{
  // unpack params
  const double*   H = (double*) grid;
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

  int i,j,idx,idx_zs,m,mm;
  double qm, h3=h*h*h;
  double sx[4] MEM_ALIGNED;
  double sy[4] MEM_ALIGNED;
  double sz[4] MEM_ALIGNED;

  // hold entire zz vector
  __m256d rZZ0, rZZ1;
  __m256d rC, rCX, rCY;
  __m256d rH0, rH1; 
  __m256d rZS0, rZS1;
  __m256d rZFZ0, rZFZ1;
  __m256d rFX, rFY, rFZ;
    
#ifdef CALC_ENERGY
  double s[4]  MEM_ALIGNED;
  __m256d rP, rCP;
#endif

  const int incrj = params->npdims[2]-8;
  const int incri = params->npdims[2]*(params->npdims[1]-8);

  for(m=0; m<N; m++){
    idx = spline->idx[m];
    mm  = spline->ind[m];
    qm  = q[mm];

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

    if(idx%4==0){ // H[idx] is 32-aligned so vectorization simple
      for(i = 0; i<8; i++){
	for(j = 0; j<8; j++){
	  rC  = _mm256_set1_pd( zx[m*8+i]*zy[m*8 + j]);
	  rCX = _mm256_set1_pd( zx[m*8+i]*zy[m*8 + j] * zfx[m*8 + i]);
	  rCY = _mm256_set1_pd( zx[m*8+i]*zy[m*8 + j] * zfy[m*8 + j]);
#ifdef CALC_ENERGY
	  rCP  = _mm256_set1_pd( zx[m*8+i]*zy[m*8 + j]);
#endif

	  rH0  = _mm256_load_pd( H+idx     );
	  rH1  = _mm256_load_pd( H+idx + 4 );
		 
	  rZS0 = _mm256_load_pd( zs + idx_zs     );
	  rZS1 = _mm256_load_pd( zs + idx_zs + 4 );
		 		    
	  rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCX),rZS0)));
	  rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCX),rZS1)));
		 

	  rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCY),rZS0)));
	  rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCY),rZS1)));
		 

	  rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
	  rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ1,rZZ1),rC),rZS1)));
		 

#ifdef CALC_ENERGY
	  rP =_mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCP),rZS0)));
	  rP =_mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCP),rZS1)));
		 
#endif

	  idx_zs +=8;
	  idx += incrj + 8;
	}
	idx += incri;
      }
    }
    else{ // H[idx] not 32-aligned, so use non-aligned loads
      for(i = 0; i<8; i++){
	for(j = 0; j<8; j++){
	  rC  = _mm256_set1_pd( zx[m*8+i]*zy[m*8 + j]);
	  rCX = _mm256_set1_pd( zx[m*8+i]*zy[m*8 + j] * zfx[m*8 + i]);
	  rCY = _mm256_set1_pd( zx[m*8+i]*zy[m*8 + j] * zfy[m*8 + j]);
#ifdef CALC_ENERGY
	  rCP = _mm256_set1_pd( zx[m*8+i]*zy[m*8 + j]);
#endif

	  rH0  = _mm256_loadu_pd( H+idx     );
	  rH1  = _mm256_loadu_pd( H+idx + 4 );
		 
	  rZS0 = _mm256_load_pd( zs + idx_zs     );
	  rZS1 = _mm256_load_pd( zs + idx_zs + 4 );
		 		    
	  rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCX),rZS0)));
	  rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCX),rZS1)));
		 
	  rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCY),rZS0)));
	  rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCY),rZS1)));
		 
	  rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
	  rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ1,rZZ1),rC),rZS1)));
		 

#ifdef CALC_ENERGY
	  rP =_mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCP),rZS0)));
	  rP =_mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCP),rZS1)));
#endif

	  idx_zs +=8;
	  idx += incrj + 8;
	}
	idx += incri;
      }
    }
    _mm256_store_pd(sx,rFX);
    _mm256_store_pd(sy,rFY);
    _mm256_store_pd(sz,rFZ);

    if(bClearF){
      force[m][XX] = 0;
      force[m][YY] = 0;
      force[m][ZZ] = 0;
    }

    force[m][XX] += -qm*scale*h3*(sx[0]+sx[1]+sx[2]+sx[3]);
    force[m][YY] += -qm*scale*h3*(sy[0]+sy[1]+sy[2]+sy[3]);
    force[m][ZZ] += -qm*scale*h3*(sz[0]+sz[1]+sz[2]+sz[3]);

#ifdef CALC_ENERGY
    _mm256_store_pd(s,rP);
    st->phi[m] = -scale*h3*(s[0]+s[1]+s[2]+s[3]);
#endif

  }
}

// -----------------------------------------------------------------------------
static void 
SE_int_split_AVX_P16_d(rvec* force, real* grid, real* q,  
		       splinedata_t *spline,
		       const SE_FGG_params* params, real scale,
		       gmx_bool bClearF)
{
  // unpack params
  const double*   H = (double*) grid;
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
  double sx[4] MEM_ALIGNED;
  double sy[4] MEM_ALIGNED;
  double sz[4] MEM_ALIGNED;

  // hold entire zz vector
  __m256d rZZ0, rZZ1, rZZ2, rZZ3;
  __m256d rC, rCX, rCY;
  __m256d rH0, rH1, rH2, rH3; 
  __m256d rZS0, rZS1, rZS2, rZS3;
  __m256d rZFZ0, rZFZ1, rZFZ2, rZFZ3;
  __m256d rFX, rFY, rFZ;
    
#ifdef CALC_ENERGY
  double s[4]  MEM_ALIGNED;
  __m256d rP, rCP;
#endif

  const int incrj = params->npdims[2]-16;
  const int incri = params->npdims[2]*(params->npdims[1]-16);

  for(m=0; m<N; m++){
    idx = spline->idx[m];
    mm  = spline->ind[m];
    qm = q[mm];

    _mm_prefetch( (void*) (H+idx), _MM_HINT_T0 );
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

    if(idx%4==0){ // H[idx] is 32-aligned so vectorization simple
      for(i = 0; i<16; i++){
	for(j = 0; j<16; j++){
	  rC  = _mm256_set1_pd( zx[m*16+i]*zy[m*16 + j]);
	  rCX = _mm256_set1_pd( zx[m*16+i]*zy[m*16 + j] * zfx[m*16 + i]);
	  rCY = _mm256_set1_pd( zx[m*16+i]*zy[m*16 + j] * zfy[m*16 + j]);
#ifdef CALC_ENERGY
	  rCP  = _mm256_set1_pd( zx[m*16+i]*zy[m*16 + j]);
#endif

	  rH0  = _mm256_load_pd( H+idx     );
	  rH1  = _mm256_load_pd( H+idx + 4 );
	  rH2  = _mm256_load_pd( H+idx + 8 );
	  rH3  = _mm256_load_pd( H+idx + 12);

	  rZS0 = _mm256_load_pd( zs + idx_zs     );
	  rZS1 = _mm256_load_pd( zs + idx_zs + 4 );
	  rZS2 = _mm256_load_pd( zs + idx_zs + 8 );
	  rZS3 = _mm256_load_pd( zs + idx_zs + 12);
		    
	  rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCX),rZS0)));
	  rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCX),rZS1)));
	  rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rCX),rZS2)));
	  rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rCX),rZS3)));

	  rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCY),rZS0)));
	  rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCY),rZS1)));
	  rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rCY),rZS2)));
	  rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rCY),rZS3)));

	  rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
	  rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ1,rZZ1),rC),rZS1)));
	  rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ2,rZZ2),rC),rZS2)));
	  rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ3,rZZ3),rC),rZS3)));

#ifdef CALC_ENERGY
	  rP =_mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCP),rZS0)));
	  rP =_mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCP),rZS1)));
	  rP =_mm256_add_pd(rP,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rCP),rZS2)));
	  rP =_mm256_add_pd(rP,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rCP),rZS3)));
#endif

	  idx_zs +=16;
	  idx += incrj + 16;
	}
	idx += incri;
      }
    }
    else{ // H[idx] not 32-aligned, so use non-aligned loads
      for(i = 0; i<16; i++){
	for(j = 0; j<16; j++){
	  rC  = _mm256_set1_pd( zx[m*16+i]*zy[m*16 + j]);
	  rCX = _mm256_set1_pd( zx[m*16+i]*zy[m*16 + j] * zfx[m*16 + i]);
	  rCY = _mm256_set1_pd( zx[m*16+i]*zy[m*16 + j] * zfy[m*16 + j]);
#ifdef CALC_ENERGY
	  rCP = _mm256_set1_pd( zx[m*16+i]*zy[m*16 + j]);
#endif

	  rH0  = _mm256_loadu_pd( H+idx     );
	  rH1  = _mm256_loadu_pd( H+idx + 4 );
	  rH2  = _mm256_loadu_pd( H+idx + 8 );
	  rH3  = _mm256_loadu_pd( H+idx + 12);

	  rZS0 = _mm256_load_pd( zs + idx_zs     );
	  rZS1 = _mm256_load_pd( zs + idx_zs + 4 );
	  rZS2 = _mm256_load_pd( zs + idx_zs + 8 );
	  rZS3 = _mm256_load_pd( zs + idx_zs + 12);
		    
	  rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCX),rZS0)));
	  rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCX),rZS1)));
	  rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rCX),rZS2)));
	  rFX =_mm256_add_pd(rFX,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rCX),rZS3)));

	  rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCY),rZS0)));
	  rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCY),rZS1)));
	  rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rCY),rZS2)));
	  rFY =_mm256_add_pd(rFY,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rCY),rZS3)));

	  rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
	  rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ1,rZZ1),rC),rZS1)));
	  rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ2,rZZ2),rC),rZS2)));
	  rFZ =_mm256_add_pd(rFZ,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rZFZ3,rZZ3),rC),rZS3)));

#ifdef CALC_ENERGY
	  rP =_mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rCP),rZS0)));
	  rP =_mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rCP),rZS1)));
	  rP =_mm256_add_pd(rP,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rCP),rZS2)));
	  rP =_mm256_add_pd(rP,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rCP),rZS3)));
#endif

	  idx_zs +=16;
	  idx += incrj + 16;
	}
	idx += incri;
      }
    }
    _mm256_store_pd(sx,rFX);
    _mm256_store_pd(sy,rFY);
    _mm256_store_pd(sz,rFZ);

    if(bClearF){
      force[m][XX] = 0;
      force[m][YY] = 0;
      force[m][ZZ] = 0;
    }

    force[m][XX] += -qm*scale*h3*(sx[0]+sx[1]+sx[2]+sx[3]);
    force[m][YY] += -qm*scale*h3*(sy[0]+sy[1]+sy[2]+sy[3]);
    force[m][ZZ] += -qm*scale*h3*(sz[0]+sz[1]+sz[2]+sz[3]);

#ifdef CALC_ENERGY
    _mm256_stream_pd(s,rP);
    st->phi[m] = -scale*h3*(s[0]+s[1]+s[2]+s[3]);
#endif

  }
}

// -----------------------------------------------------------------------------
static void 
SE_int_split_AVX_dispatch_d(rvec* force, real* grid, real* q, 
			    splinedata_t *spline,
			    const SE_FGG_params* params, real scale,
			    gmx_bool bClearF)
{
  const int p = params->P;
  const int incrj = params->dims[2]; // middle increment
  const int incri = params->npdims[2]*(params->dims[1]);// outer increment
  
  // if P, incri or increments are not divisible by 4, fall back on vanilla
  if( isnot_div_by_4(p) || isnot_div_by_4(incri) || isnot_div_by_4(incrj) ){
    __DISPATCHER_MSG("[FGG INT AVX DOUBLE] AVX Abort (PARAMS)\n");
    SE_int_split_SSE_dispatch_d(force, grid, q, spline, params, scale, bClearF);
    return;
  }
    
  // otherwise the preconditions for AVX codes are satisfied. 
  if(p==8){
    // specific for p=8
    __DISPATCHER_MSG("[FGG INT AVX DOUBLE] P=8\n");
    SE_int_split_AVX_P8_d(force, grid, q, spline, params, scale, bClearF);
  }
  else if(p==16){
    // specific for p=16
    __DISPATCHER_MSG("[FGG INT AVX DOUBLE] P=16\n");
    SE_int_split_AVX_P16_d(force, grid, q, spline, params, scale, bClearF);
  }
  else if(p%8==0){
    // for p divisible by 8
    __DISPATCHER_MSG("[FGG INT AVX DOUBLE] P unroll 8\n");
    SE_int_split_AVX_u8_d(force, grid, q, spline, params, scale, bClearF); 
  }
  else if(p%4==0){
    // vanilla AVX code (p divisible by 4)
    __DISPATCHER_MSG("[FGG INT AVX DOUBLE] P unroll 4\n");
    SE_int_split_AVX_d(force, grid, q, spline, params, scale, bClearF);
  }
  else {
    // vanilla SSE code (any even p)
    __DISPATCHER_MSG("[FGG INT AVX DOUBLE] Vanilla\n");
    SE_int_split_SSE_d(force, grid, q, spline, params, scale, bClearF);
  }
}

#endif //GMX_SIMD_X86_AVX_256
#endif //GMX_DOUBLE
#endif //_SE_INT_AVX_256_DOUBLE_
