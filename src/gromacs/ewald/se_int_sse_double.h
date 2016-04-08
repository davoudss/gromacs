#ifndef _SE_INT_SSE_DOUBLE_
#define _SE_INT_SSE_DOUBLE_

/* SE int SSE double integration */
#if GMX_DOUBLE==1
#include "se.h"

static void 
SE_int_split_d(rvec *force, real *grid, real *q,
	       splinedata_t *spline,
	       const SE_FGG_params *params, real scale,
	       gmx_bool bClearF) 
{
  // unpack parameters
  const real*    H = (real*) grid;
  const real*   zs = (real*) spline->zs;
  const real*   zx = (real*) spline->theta[0];
  const real*   zy = (real*) spline->theta[1];
  const real*   zz = (real*) spline->theta[2];
  const real*   zfx= (real*) spline->dtheta[0];
  const real*   zfy= (real*) spline->dtheta[1];
  const real*   zfz= (real*) spline->dtheta[2];

  const int    p = params->P;
  const double h = params->h;

  int i, j, k, m, mm, idx, idx_zs, idx_zz;
  real force_m[3], cij, Hzc, qm, h3=h*h*h;
#ifdef CALC_ENERGY
  real phi_m;
#endif

  const int incrj = params->npdims[2]-p; // middle increment
  const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

#ifdef _OPENMP
#pragma omp for private(mm) schedule(static) // work-share over OpenMP threads here
#endif
  for(mm=0; mm<spline->n; mm++)
    {
      idx = spline->idx[mm];
      m   = spline->ind[mm];
      qm  = q[m];
      force_m[0] = 0; force_m[1] = 0; force_m[2] = 0;
#ifdef CALC_ENERGY
      phi_m = 0;
#endif

      idx_zs = 0;
      for(i = 0; i<p; i++)
        {
          for(j = 0; j<p; j++)
            {
              cij = zx[p*m+i]*zy[p*m+j];
              idx_zz=p*m;
              for(k = 0; k<p; k++)
                {
		  Hzc    = H[idx]*zs[idx_zs]*zz[idx_zz]*cij;
#ifdef CALC_ENERGY
		  phi_m += Hzc;
#endif
		  force_m[0] += Hzc*zfx[m*p+i];
		  force_m[1] += Hzc*zfy[m*p+j];
		  force_m[2] += Hzc*zfz[m*p+k];
		  
                  idx++; idx_zs++; idx_zz++;
                }
              idx += incrj;
            }
          idx += incri;
        }
      if(bClearF)
	{
	  force[m][XX] = 0;
	  force[m][YY] = 0;
	  force[m][ZZ] = 0;
	}

      force[m][XX] += -qm*h3*scale*force_m[0];
      force[m][YY] += -qm*h3*scale*force_m[1];
      force[m][ZZ] += -qm*h3*scale*force_m[2];

#ifdef CALC_ENERGY
      st->phi[m]   = -h3*scale*phi_m;
#endif
    }
}


// -----------------------------------------------------------------------------
static void 
SE_int_split_SSE_d(rvec *force, real *grid, real *q,
		   splinedata_t *spline,
		   const SE_FGG_params *params, real scale,
		   gmx_bool bClearF)
{
  // unpack params
  const double*    H = (double*) grid;
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

  int i,j,k,m,idx,idx_zs,idx_zz, mm,mp,imp,jmp,kmp;
  double qm, h3 = h*h*h;
  double sx[2] MEM_ALIGNED;
  double sy[2] MEM_ALIGNED;
  double sz[2] MEM_ALIGNED;
  
  __m128d rH0, rZZ0, rZS0,rZFZ0;
  __m128d rC, rCX, rCY;
  __m128d rFX, rFY, rFZ;
#ifdef CALC_ENERGY
  double s[2]  MEM_ALIGNED;
  __m128d rP, rCP;
#endif

  const int incrj = params->npdims[2]-p;
  const int incri = params->npdims[2]*(params->npdims[1]-p);

  for(mm=0; mm<N; mm++){
    idx = spline->idx[mm];
    m = spline->ind[mm];
    qm = q[m];
    mp = m*p;
    
    idx_zs = 0;
    rFX = _mm_setzero_pd();
    rFY = _mm_setzero_pd();
    rFZ = _mm_setzero_pd();
    
#ifdef CALC_ENERGY
    rP = _mm_setzero_pd();
#endif
    
    if(idx%2==0){ // H[idx] is 16-aligned so vectorization simple
      for(i = 0; i<p; i++){
	imp = mp+i;
	for(j = 0; j<p; j++){
	  jmp = mp+j;
	  rC  = _mm_set1_pd( zx[imp]*zy[jmp] );
	  rCX = _mm_set1_pd(zx[imp]*zy[jmp]*zfx[imp]);
	  rCY = _mm_set1_pd(zx[imp]*zy[jmp]*zfy[jmp]);
#ifdef CALC_ENERGY
	  rCP  = _mm_set1_pd( zx[imp]*zy[jmp]);
#endif
	  
	  idx_zz=mp;
	  for(k = 0; k<p; k+=2){
	    kmp = mp+k;
	    rZFZ0= _mm_load_pd( zfz+ kmp );
	    rH0  = _mm_load_pd( H  + idx );
	    rZZ0 = _mm_load_pd( zz + idx_zz);
	    rZS0 = _mm_load_pd( zs + idx_zs);
	    
	    rH0 = _mm_mul_pd(rH0,_mm_mul_pd(rZZ0,rZS0));
	    rFX = _mm_add_pd(rFX,_mm_mul_pd(rH0,rCX));
	    rFY = _mm_add_pd(rFY,_mm_mul_pd(rH0,rCY));
	    rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(rC,rZFZ0)));
	    
#ifdef CALC_ENERGY
	    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,rCP));
#endif			
	    idx+=2; 
	    idx_zs+=2; 
	    idx_zz+=2;
	  }
	  idx += incrj;
	}
	idx += incri;
      }
    }
    else { // H[idx] not 16-aligned, so use non-aligned loads
      for(i = 0; i<p; i++){
	imp = mp+i;
	for(j = 0; j<p; j++){
	  jmp = mp+j;
	  rC  = _mm_set1_pd( zx[imp]*zy[jmp] );
	  rCX = _mm_set1_pd(zx[imp]*zy[jmp]*zfx[imp]);
	  rCY = _mm_set1_pd(zx[imp]*zy[jmp]*zfy[jmp]);
#ifdef CALC_ENERGY
	  rCP = _mm_set1_pd( zx[imp]*zy[jmp] );
#endif
	  idx_zz=m*p;
	  for(k = 0; k<p; k+=2){
	    kmp = mp+k;
	    rZFZ0= _mm_load_pd( zfz + kmp );
	    rH0  = _mm_loadu_pd( H+idx );
	    rZZ0 = _mm_load_pd( zz + idx_zz);
	    rZS0 = _mm_load_pd( zs + idx_zs);
	    
	    rH0 = _mm_mul_pd(rH0,_mm_mul_pd(rZZ0,rZS0));
	    rFX = _mm_add_pd(rFX,_mm_mul_pd(rH0,rCX));
	    rFY = _mm_add_pd(rFY,_mm_mul_pd(rH0,rCY));
	    rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(rC,rZFZ0)));
	    
#ifdef CALC_ENERGY
	    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,rCP));
#endif			
	    idx+=2; 
	    idx_zs+=2; 
	    idx_zz+=2;
	  }
	  idx += incrj;
	}
	idx += incri;
      }
      
    }
    _mm_store_pd(sx,rFX);
    _mm_store_pd(sy,rFY);
    _mm_store_pd(sz,rFZ);
    
    force[m][XX] = -qm*scale*h3*(sx[0]+sx[1]);
    force[m][YY] = -qm*scale*h3*(sy[0]+sy[1]);
    force[m][ZZ] = -qm*scale*h3*(sz[0]+sz[1]);
    
#ifdef CALC_ENERGY
    _mm_store_pd(s,rP);
    st->phi[m] = -scale*h3*(s[0]+s[1]);
#endif
  }
}


// -----------------------------------------------------------------------------
static void 
SE_int_split_SSE_u8_d(rvec *force, real *grid, real *q,
		      splinedata_t *spline,
		      const SE_FGG_params *params, real scale,
		      gmx_bool bClearF)
{
  // unpack params
  const double*    H = (double*) grid;
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
  double sx[2] MEM_ALIGNED;
  double sy[2] MEM_ALIGNED;
  double sz[2] MEM_ALIGNED;

  __m128d rH0, rZZ0, rZS0, rZFZ0;
  __m128d rH1, rZZ1, rZS1, rZFZ1;
  __m128d rH2, rZZ2, rZS2, rZFZ2;
  __m128d rH3, rZZ3, rZS3, rZFZ3;
  __m128d rFX, rFY, rFZ;
  __m128d  rC, rCX, rCY;

#ifdef CALC_ENERGY
  double s[2]  MEM_ALIGNED;
  __m128d rP, rCP;
#endif

  const int incrj = params->npdims[2]-p;
  const int incri = params->npdims[2]*(params->npdims[1]-p);

  for(m=0; m<N; m++){
    idx = spline->idx[m];
    mm  = spline->ind[m];
    qm = q[mm];
    
    _mm_prefetch( (void*) (H+idx), _MM_HINT_T0);
    
    idx_zs = 0;
    _mm_prefetch( (void*) zs, _MM_HINT_T0);
    
    rFX = _mm_setzero_pd();
    rFY = _mm_setzero_pd();
    rFZ = _mm_setzero_pd();
#ifdef CALC_ENERGY
    rP  = _mm_setzero_pd();
#endif
    
    if(idx%2==0){ // H[idx] is 16-aligned so vectorization simple
      for(i = 0; i<p; i++){
	for(j = 0; j<p; j++){
	  rC  = _mm_set1_pd( zx[m*p+i]*zy[m*p+j] );
	  rCX = _mm_set1_pd(zx[m*p+i]*zy[m*p+j]*zfx[m*p+i]);
	  rCY = _mm_set1_pd(zx[m*p+i]*zy[m*p+j]*zfy[m*p+j]);
#ifdef CALC_ENERGY
	  rCP = _mm_set1_pd( zx[m*p+i]*zy[m*p+j] );
#endif
	  
	  idx_zz=m*p;
	  
	  for(k = 0; k<p; k+=8){
	    rH0  = _mm_load_pd( H+idx    );
	    rH1  = _mm_load_pd( H+idx + 2);
	    rH2  = _mm_load_pd( H+idx + 4);
	    rH3  = _mm_load_pd( H+idx + 6);
	    
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
	    
	    rFX = _mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCX),rZS0)));
	    rFX = _mm_add_pd(rFX,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rCX),rZS1)));
	    rFX = _mm_add_pd(rFX,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rCX),rZS2)));
	    rFX = _mm_add_pd(rFX,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rCX),rZS3)));
	    
	    rFY = _mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCY),rZS0)));
	    rFY = _mm_add_pd(rFY,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rCY),rZS1)));
	    rFY = _mm_add_pd(rFY,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rCY),rZS2)));
	    rFY = _mm_add_pd(rFY,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rCY),rZS3)));
	    
	    rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
	    rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ1,rZZ1),rC),rZS1)));
	    rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ2,rZZ2),rC),rZS2)));
	    rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ3,rZZ3),rC),rZS3)));

#ifdef CALC_ENERGY
	    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCP),rZS0)));
	    rP = _mm_add_pd(rP,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rCP),rZS1)));
	    rP = _mm_add_pd(rP,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rCP),rZS2)));
	    rP = _mm_add_pd(rP,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rCP),rZS3)));
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
    else{ // H[idx] not 16-aligned, so use non-aligned load from H
      for(i = 0; i<p; i++){
	for(j = 0; j<p; j++){
	  rC  = _mm_set1_pd( zx[m*p+i]*zy[m*p+j] );
	  rCX = _mm_set1_pd(zx[m*p+i]*zy[m*p+j]*zfx[m*p+i]);
	  rCY = _mm_set1_pd(zx[m*p+i]*zy[m*p+j]*zfy[m*p+j]);
#ifdef CALC_ENERGY
	  rCP = _mm_set1_pd( zx[m*p+i]*zy[m*p+j] );
#endif
	  
	  idx_zz=m*p;
	  for(k = 0; k<p; k+=8){
	    rH0  = _mm_loadu_pd( H+idx    );
	    rH1  = _mm_loadu_pd( H+idx + 2);
	    rH2  = _mm_loadu_pd( H+idx + 4);
	    rH3  = _mm_loadu_pd( H+idx + 6);
	    
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

	    rFX = _mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCX),rZS0)));
	    rFX = _mm_add_pd(rFX,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rCX),rZS1)));
	    rFX = _mm_add_pd(rFX,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rCX),rZS2)));
	    rFX = _mm_add_pd(rFX,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rCX),rZS3)));

	    rFY = _mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCY),rZS0)));
	    rFY = _mm_add_pd(rFY,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rCY),rZS1)));
	    rFY = _mm_add_pd(rFY,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rCY),rZS2)));
	    rFY = _mm_add_pd(rFY,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rCY),rZS3)));
			
	    rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
	    rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ1,rZZ1),rC),rZS1)));
	    rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ2,rZZ2),rC),rZS2)));
	    rFZ = _mm_add_pd(rFZ,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ3,rZZ3),rC),rZS3)));

#ifdef CALC_ENERGY
	    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCP),rZS0)));
	    rP = _mm_add_pd(rP,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rCP),rZS1)));
	    rP = _mm_add_pd(rP,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rCP),rZS2)));
	    rP = _mm_add_pd(rP,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rCP),rZS3)));
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
    _mm_store_pd(sx,rFX);
    _mm_store_pd(sy,rFY);
    _mm_store_pd(sz,rFZ);
    
    force[m][XX] = -qm*scale*h3*(sx[0]+sx[1]);
    force[m][YY] = -qm*scale*h3*(sy[0]+sy[1]);
    force[m][ZZ] = -qm*scale*h3*(sz[0]+sz[1]);
    
#ifdef CALC_ENERGY
    _mm_store_pd(s,rP);
    st->phi[m] = -scale*h3*(s[0]+s[1]);
#endif
  }
}

// -----------------------------------------------------------------------------
static void 
SE_int_split_SSE_P8_d(rvec *force, real *grid, real *q,
		      splinedata_t *spline,
		      const SE_FGG_params *params, real scale,
		      gmx_bool bClearF)
{
  // unpack params
  const double*    H = (double*) grid;
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

  int i,j, m, idx,idx_zs, mm,m8,im8,jm8;
  double qm, h3=h*h*h;
  double sx[2] MEM_ALIGNED;
  double sy[2] MEM_ALIGNED;
  double sz[2] MEM_ALIGNED;

  // hold entire zz vector
  __m128d rZZ0, rZZ1, rZZ2, rZZ3;
  __m128d rC, rCX, rCY;
  __m128d rH0, rH1, rH2, rH3; 
  __m128d rZS0, rZS1, rZS2, rZS3;
  __m128d rZFZ0, rZFZ1, rZFZ2, rZFZ3;
  __m128d rFX, rFY, rFZ;

#ifdef CALC_ENERGY
  double s[2]  MEM_ALIGNED;
  __m128d rP, rCP;
#endif

  const int incrj = params->npdims[2]-8;
  const int incri = params->npdims[2]*(params->npdims[1]-8);

  for(mm=0; mm<N; mm++){
    idx = spline->idx[mm];
    m = spline->ind[mm];
    qm = q[m];
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
    
    if(idx%2==0){ // H[idx] is 16-aligned so vectorization simple
      for(i = 0; i<8; i++){
	im8 = i + m8;
	for(j = 0; j<8; j++){
	  jm8 = j + m8;
	  rC  = _mm_set1_pd( zx[im8]*zy[jm8]);
	  rCX = _mm_set1_pd( zx[im8]*zy[jm8] * zfx[im8]);
	  rCY = _mm_set1_pd( zx[im8]*zy[jm8] * zfy[jm8]);
#ifdef CALC_ENERGY
	  rCP  = _mm_set1_pd( zx[im8]*zy[jm8]);
#endif
	  
	  rH0  = _mm_load_pd( H+idx    );
	  rH1  = _mm_load_pd( H+idx + 2);
	  rH2  = _mm_load_pd( H+idx + 4);
	  rH3  = _mm_load_pd( H+idx + 6);
		  
	  rZS0 = _mm_load_pd( zs + idx_zs    );
	  rZS1 = _mm_load_pd( zs + idx_zs + 2);
	  rZS2 = _mm_load_pd( zs + idx_zs + 4);
	  rZS3 = _mm_load_pd( zs + idx_zs + 6);

	  rH0 = _mm_mul_pd(rH0,_mm_mul_pd(rZZ0,rZS0));
	  rH1 = _mm_mul_pd(rH1,_mm_mul_pd(rZZ1,rZS1));
	  rH2 = _mm_mul_pd(rH2,_mm_mul_pd(rZZ2,rZS2));
	  rH3 = _mm_mul_pd(rH3,_mm_mul_pd(rZZ3,rZS3));

	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,rCX));
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH1,rCX));
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH2,rCX));
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH3,rCX));

	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,rCY));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH1,rCY));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH2,rCY));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH3,rCY));

	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(rZFZ0,rC)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH1,_mm_mul_pd(rZFZ1,rC)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH2,_mm_mul_pd(rZFZ2,rC)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH3,_mm_mul_pd(rZFZ3,rC)));

#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,rCP));
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH1,rCP));
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH2,rCP));
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH3,rCP));
#endif

	  idx_zs +=8;
	  idx += incrj + 8;
	}
	idx += incri;
      }
    }
    else{ // H[idx] not 16-aligned, so use non-aligned loads
      for(i = 0; i<8; i++){
	im8 = i + m8;
	for(j = 0; j<8; j++){
	  jm8 = j + m8;
	  rC  = _mm_set1_pd( zx[im8]*zy[jm8]);
	  rCX = _mm_set1_pd( zx[im8]*zy[jm8] * zfx[im8]);
	  rCY = _mm_set1_pd( zx[im8]*zy[jm8] * zfy[jm8]);
#ifdef CALC_ENERGY
	  rCP = _mm_set1_pd( zx[im8]*zy[jm8]);
#endif

	  rH0  = _mm_loadu_pd( H+idx    );
	  rH1  = _mm_loadu_pd( H+idx + 2);
	  rH2  = _mm_loadu_pd( H+idx + 4);
	  rH3  = _mm_loadu_pd( H+idx + 6);

	  rZS0 = _mm_load_pd( zs + idx_zs    );
	  rZS1 = _mm_load_pd( zs + idx_zs + 2);
	  rZS2 = _mm_load_pd( zs + idx_zs + 4);
	  rZS3 = _mm_load_pd( zs + idx_zs + 6);

	  rH0 = _mm_mul_pd(rH0,_mm_mul_pd(rZZ0,rZS0));
	  rH1 = _mm_mul_pd(rH1,_mm_mul_pd(rZZ1,rZS1));
	  rH2 = _mm_mul_pd(rH2,_mm_mul_pd(rZZ2,rZS2));
	  rH3 = _mm_mul_pd(rH3,_mm_mul_pd(rZZ3,rZS3));

	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,rCX));
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH1,rCX));
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH2,rCX));
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH3,rCX));

	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,rCY));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH1,rCY));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH2,rCY));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH3,rCY));

	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(rZFZ0,rC)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH1,_mm_mul_pd(rZFZ1,rC)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH2,_mm_mul_pd(rZFZ2,rC)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH3,_mm_mul_pd(rZFZ3,rC)));

#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,rCP));
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH1,rCP));
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH2,rCP));
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH3,rCP));
#endif

	  idx_zs +=8;
	  idx += incrj + 8;
	}
	idx += incri;
      }
    }
    _mm_store_pd(sx,rFX);
    _mm_store_pd(sy,rFY);
    _mm_store_pd(sz,rFZ);
      
    force[m][XX] = -qm*scale*h3*(sx[0]+sx[1]);
    force[m][YY] = -qm*scale*h3*(sy[0]+sy[1]);
    force[m][ZZ] = -qm*scale*h3*(sz[0]+sz[1]);
      
#ifdef CALC_ENERGY
    _mm_store_pd(s,rP);
    st->phi[m] = -scale*h3*(s[0]+s[1]);
#endif
  }
}

// -----------------------------------------------------------------------------
static void 
SE_int_split_SSE_P16_d(rvec *force, real *grid, real *q,
		       splinedata_t *spline,
		       const SE_FGG_params *params, real scale,
		       gmx_bool bClearF)
{
  // unpack params
  const double*    H = (double*) grid;
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
  double sx[2] MEM_ALIGNED;
  double sy[2] MEM_ALIGNED;
  double sz[2] MEM_ALIGNED;
    

  // hold entire zz vector
  __m128d rZZ0 , rZZ1 , rZZ2 , rZZ3 , rZZ4 , rZZ5 , rZZ6 , rZZ7; 
  __m128d rZFZ0, rZFZ1, rZFZ2, rZFZ3, rZFZ4, rZFZ5, rZFZ6, rZFZ7;
  __m128d rC, rCX, rCY, rFX, rFY, rFZ;
  __m128d rH0, rZS0;

#ifdef CALC_ENERGY
  double s[2]  MEM_ALIGNED;
  __m128d rP, rCP;
#endif

  const int incrj = params->npdims[2]-16;
  const int incri = params->npdims[2]*(params->npdims[1]-16);

  for(mm=0; mm<N; mm++){
    idx = spline->idx[mm];
    m = spline->ind[mm];
    qm = q[m];
    m16 = m*16;
    
    _mm_prefetch( (void*) (H+idx), _MM_HINT_T0);
    
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

    if(idx%2==0){ // H[idx] is 16-aligned so vectorization simple
      for(i = 0; i<16; i++){
	im16 = i + m16;
	for(j = 0; j<16; j++){
	  jm16 = j + m16;
	  rC  = _mm_set1_pd( zx[im16]*zy[jm16]);
	  rCX = _mm_set1_pd( zx[im16]*zy[jm16]*zfx[im16]);
	  rCY = _mm_set1_pd( zx[im16]*zy[jm16]*zfy[jm16]);
#ifdef CALC_ENERGY
	  rCP = _mm_set1_pd( zx[im16]*zy[jm16]);
#endif
	  
	  /* 0 */ 
	  rH0  = _mm_load_pd( H+idx );
	  rZS0 = _mm_load_pd( zs + idx_zs);		    
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
	  
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ0,rZZ0),rCP),rZS0)));
#endif
	  /* 1 */ 
	  rH0  = _mm_load_pd( H+idx + 2);
	  rZS0 = _mm_load_pd( zs + idx_zs + 2);		    
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ1,rZZ1),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rCP),rZS0)));
#endif
	  
	  /* 2 */ 
	  rH0  = _mm_load_pd( H+idx + 4);
	  rZS0 = _mm_load_pd( zs + idx_zs + 4);		    
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ2,rZZ2),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rCP),rZS0)));
#endif
	  
	  /* 3 */ 
	  rH0  = _mm_load_pd( H+idx + 6);
	  rZS0 = _mm_load_pd( zs + idx_zs + 6);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ3,rZZ3),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rCP),rZS0)));
#endif

	  /* 4 */ 
	  rH0  = _mm_load_pd( H+idx + 8);
	  rZS0 = _mm_load_pd( zs + idx_zs + 8);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ4,rZZ4),rC),rZS0)));
#ifdef CALC_NERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rCP),rZS0)));
#endif

	  /* 5 */ 
	  rH0  = _mm_load_pd( H+idx + 10);
	  rZS0 = _mm_load_pd( zs + idx_zs + 10);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ5,rZZ5),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rCP),rZS0)));
#endif

	  /* 6 */ 
	  rH0  = _mm_load_pd( H+idx + 12);
	  rZS0 = _mm_load_pd( zs + idx_zs + 12);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ6,rZZ6),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rCP),rZS0)));
#endif
	  /* 7 */ 
	  rH0  = _mm_load_pd( H+idx + 14);
	  rZS0 = _mm_load_pd( zs + idx_zs + 14);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ7,rZZ7),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rCP),rZS0)));
#endif

	  idx_zs +=16;
	  idx += incrj + 16;
	}
	idx += incri;
      }
    }
    else{ // H[idx] not 16-aligned, so use non-aligned loads
      for(i = 0; i<16; i++){
	im16 = i + m16;
	for(j = 0; j<16; j++){
	  jm16 = j + m16;
	  rC  = _mm_set1_pd( zx[im16]*zy[jm16]);
	  rCX = _mm_set1_pd( zx[im16]*zy[jm16]*zfx[im16]);
	  rCY = _mm_set1_pd( zx[im16]*zy[jm16]*zfy[jm16]);
#ifdef CALC_ENERGY
	  rCP  = _mm_set1_pd( zx[im16]*zy[jm16]);
#endif

	  /* 0 */ 
	  rH0  = _mm_loadu_pd( H+idx );
	  rZS0 = _mm_load_pd( zs + idx_zs);		    
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ0,rZZ0),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rCP),rZS0)));
#endif

	  /* 1 */ 
	  rH0  = _mm_loadu_pd( H+idx + 2);
	  rZS0 = _mm_load_pd( zs + idx_zs + 2);		    
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ1,rZZ1),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rCP),rZS0)));
#endif

	  /* 2 */ 
	  rH0  = _mm_loadu_pd( H+idx + 4);
	  rZS0 = _mm_load_pd( zs + idx_zs + 4);		    
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ2,rZZ2),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rCP),rZS0)));
#endif

	  /* 3 */ 
	  rH0  = _mm_loadu_pd( H+idx + 6);
	  rZS0 = _mm_load_pd( zs + idx_zs + 6);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ3,rZZ3),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rCP),rZS0)));
#endif

	  /* 4 */ 
	  rH0  = _mm_loadu_pd( H+idx + 8);
	  rZS0 = _mm_load_pd( zs + idx_zs + 8);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ4,rZZ4),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rCP),rZS0)));
#endif

	  /* 5 */ 
	  rH0  = _mm_loadu_pd( H+idx + 10);
	  rZS0 = _mm_load_pd( zs + idx_zs + 10);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ5,rZZ5),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rCP),rZS0)));
#endif

	  /* 6 */ 
	  rH0  = _mm_loadu_pd( H+idx + 12);
	  rZS0 = _mm_load_pd( zs + idx_zs + 12);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ6,rZZ6),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rCP),rZS0)));
#endif

	  /* 7 */ 
	  rH0  = _mm_loadu_pd( H+idx + 14);
	  rZS0 = _mm_load_pd( zs + idx_zs + 14);
	  rFX =_mm_add_pd(rFX,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rCX),rZS0)));
	  rFY =_mm_add_pd(rFY,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rCY),rZS0)));
	  rFZ =_mm_add_pd(rFZ,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(_mm_mul_pd(rZFZ7,rZZ7),rC),rZS0)));
#ifdef CALC_ENERGY
	  rP =_mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rCP),rZS0)));
#endif
	  idx_zs +=16;
	  idx += incrj + 16;
	}
	idx += incri;
      }
    }

    _mm_store_pd(sx,rFX);
    _mm_store_pd(sy,rFY);
    _mm_store_pd(sz,rFZ);

    force[m][XX] = -qm*scale*h3*(sx[0]+sx[1]);
    force[m][YY] = -qm*scale*h3*(sy[0]+sy[1]);
    force[m][ZZ] = -qm*scale*h3*(sz[0]+sz[1]);

#ifdef CALC_ENERGY
    _mm_store_pd(s,rP);
    st->phi[m] = -scale*h3*(s[0]+s[1]);
#endif
  }
}

// -----------------------------------------------------------------------------
static void 
SE_int_split_SSE_dispatch_d(rvec *force, real *grid, real *q,
			    splinedata_t *spline,
			    const SE_FGG_params *params, real scale,
			    gmx_bool bClearF)
{
  const int p = params->P;
  const int incrj = params->dims[2]; // middle increment
  const int incri = params->npdims[2]*(params->dims[1]);// outer increment

  // if P is odd, or if either increment is odd, fall back on vanilla
  if( is_odd(p) || is_odd(incri) || is_odd(incrj) ){
    __DISPATCHER_MSG("[FGG INT SSE DOUBLE] SSE Abort (PARAMS)\n");
    SE_int_split_d(force, grid, q, spline, params, scale, bClearF);
    return;
  }
  // otherwise the preconditions for SSE codes are satisfied. 
  if(p==8){
    // specific for p=8
    __DISPATCHER_MSG("[FGG INT SSE DOUBLE] P=8\n");
    SE_int_split_SSE_P8_d(force, grid, q, spline, params, scale, bClearF);
  }
  else if(p==16){
    // specific for p=16
    __DISPATCHER_MSG("[FGG INT SSE DOUBLE] P=16\n");
    SE_int_split_SSE_P16_d(force, grid, q, spline, params, scale, bClearF); 
  }
  else if(p%8==0){
    // for p divisible by 8
    __DISPATCHER_MSG("[FGG INT SSE DOUBLE] P unroll 8\n");
    SE_int_split_SSE_u8_d(force, grid, q, spline, params, scale, bClearF); 
  }
  else{
    // vanilla SSE code (any even p)
    __DISPATCHER_MSG("[FGG INT SSE DOUBLE] Vanilla\n");
    SE_int_split_SSE_d(force, grid, q, spline, params, scale, bClearF);
  }
}


#endif // GMX_DOUBLE
#endif // _SE_INT_SSE_DOUBLE_
