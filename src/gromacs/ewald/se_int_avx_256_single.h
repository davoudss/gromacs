#ifndef _SE_INT_AVX_256_SINGLE_
#define _SE_INT_AVX_256_SINGLE_


/* SE AVX 256 single integration */
#ifdef GMX_X86_AVX_256
#include "se.h"

// -----------------------------------------------------------------------------
static void SE_int_split_AVX(rvec* force,  real* grid, real* q,
			     splinedata_t *spline,
			     const SE_FGG_params* params, real scale,
			     gmx_bool bClearF)
{
  // unpack params
  const float*   H = (float*) grid;
  const float*  zs = (float*) spline->zs;
  const float*  zx = (float*) spline->theta[0];
  const float*  zy = (float*) spline->theta[1];
  const float*  zz = (float*) spline->theta[2];
  const float* zfx = (float*) spline->dtheta[0];
  const float* zfy = (float*) spline->dtheta[1];
  const float* zfz = (float*) spline->dtheta[2];

  /* ASSUME P=8 const int p = params->P; */
  const int N = params->N;
  const float h=params->h;

  int i,j,idx,idx_zs,m,mm;
  float qm, h3=h*h*h;
  float sx[8] MEM_ALIGNED;
  float sy[8] MEM_ALIGNED;
  float sz[8] MEM_ALIGNED;

  // hold entire zz vector
  __m256 rZZ0;
  __m256 rC, rCX, rCY;
  __m256 rH0;
  __m256 rZS0;
  __m256 rZFZ0;
  __m256 rFX, rFY, rFZ;
    
#ifdef CALC_ENERGY
  float s[8]  MEM_ALIGNED;
  __m256 rP, rCP;
#endif

  const int incrj = params->npdims[2]-8;
  const int incri = params->npdims[2]*(params->npdims[1]-8);

  for(m=0; m<N; m++){
    idx = spline->idx[m];
    mm  = spline->ind[m];
    qm = q[mm];

    idx_zs = 0;
    rFX = _mm256_setzero_ps();
    rFY = _mm256_setzero_ps();
    rFZ = _mm256_setzero_ps();
#ifdef CALC_ENERGY
    rP  = _mm256_setzero_ps();
#endif

    /* hoist load of ZZ vector */
    rZZ0 = _mm256_load_ps(zz + m*8     );

    /* hoist load of ZFZ vector */
    rZFZ0 = _mm256_load_ps(zfz + m*8     );

    if(idx%8==0){ // H[idx] is 32-aligned so vectorization simple
      for(i = 0; i<8; i++){
	for(j = 0; j<8; j++){
	  rC  = _mm256_set1_ps( zx[m*8+i]*zy[m*8 + j]*qm);
	  rCX = _mm256_set1_ps( zx[m*8+i]*zy[m*8 + j] * zfx[m*8 + i]*qm);
	  rCY = _mm256_set1_ps( zx[m*8+i]*zy[m*8 + j] * zfy[m*8 + j]*qm);
#ifdef CALC_ENERGY
	  rCP  = _mm256_set1_ps( zx[m*8+i]*zy[m*8 + j]);
#endif

	  rH0  = _mm256_load_ps( H+idx     );
	  rZS0 = _mm256_load_ps( zs + idx_zs     );

	  rFX =_mm256_add_ps(rFX,_mm256_mul_ps(rH0,_mm256_mul_ps(_mm256_mul_ps(rZZ0,rCX),rZS0))); 
	  rFY =_mm256_add_ps(rFY,_mm256_mul_ps(rH0,_mm256_mul_ps(_mm256_mul_ps(rZZ0,rCY),rZS0)));
	  rFZ =_mm256_add_ps(rFZ,_mm256_mul_ps(rH0,_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(rZFZ0,rZZ0),rC),rZS0))); 

#ifdef CALC_ENERGY
	  rP =_mm256_add_ps(rP,_mm256_mul_ps(rH0,_mm256_mul_ps(_mm256_mul_ps(rZZ0,rCP),rZS0)));
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
	  rC  = _mm256_set1_ps( zx[m*8+i]*zy[m*8 + j]*qm);
	  rCX = _mm256_set1_ps( zx[m*8+i]*zy[m*8 + j] * zfx[m*8 + i]*qm);
	  rCY = _mm256_set1_ps( zx[m*8+i]*zy[m*8 + j] * zfy[m*8 + j]*qm);
#ifdef CALC_ENERGY
	  rCP = _mm256_set1_ps( zx[m*8+i]*zy[m*8 + j]);
#endif

	  rH0  = _mm256_loadu_ps( H+idx     );
	  rZS0 = _mm256_load_ps( zs + idx_zs);
		 		    
	  rFX =_mm256_add_ps(rFX,_mm256_mul_ps(rH0,_mm256_mul_ps(_mm256_mul_ps(rZZ0,rCX),rZS0)));
	  rFY =_mm256_add_ps(rFY,_mm256_mul_ps(rH0,_mm256_mul_ps(_mm256_mul_ps(rZZ0,rCY),rZS0)));
	  rFZ =_mm256_add_ps(rFZ,_mm256_mul_ps(rH0,_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(rZFZ0,rZZ0),rC),rZS0)));
	 
#ifdef CALC_ENERGY
	  rP =_mm256_add_ps(rP,_mm256_mul_ps(rH0,_mm256_mul_ps(_mm256_mul_ps(rZZ0,rCP),rZS0)));
#endif

	  idx_zs +=8;
	  idx += incrj + 8;
	}
	idx += incri;
      }
    }
    _mm256_store_ps(sx,rFX);
    _mm256_store_ps(sy,rFY);
    _mm256_store_ps(sz,rFZ);

    force[m][XX] = -scale*h3*(sx[0]+sx[1]+sx[2]+sx[3]+sx[4]+sx[5]+sx[6]+sx[7]);
    force[m][YY] = -scale*h3*(sy[0]+sy[1]+sy[2]+sy[3]+sy[4]+sy[5]+sy[6]+sy[7]);
    force[m][ZZ] = -scale*h3*(sz[0]+sz[1]+sz[2]+sz[3]+sz[4]+sz[5]+sz[6]+sz[7]);

#ifdef CALC_ENERGY
    _mm256_store_ps(s,rP);
    st->phi[m] = -scale*h3*(s[0]+s[1]+s[2]+s[3]+s[4]+s[5]+s[6]+s[7]);
#endif	
  }
}

// ----------------------------------------------------------------------------
static void 
SE_int_split_AVX_dispatch(rvec* force, real* grid, real* q,
			  splinedata_t *spline,
			  const SE_FGG_params* params, real scale,
			  gmx_bool bClearF)
{
  const int p = params->P;
  const int incrj = params->dims[2]; // middle increment
  const int incri = params->npdims[2]*(params->dims[1]);// outer increment

  // if P, incri or increments are not divisible by 4, fall back on vanilla
  if( isnot_div_by_8(p) || isnot_div_by_8(incri) || isnot_div_by_8(incrj) || (p%8)!=0)
    {
      __DISPATCHER_MSG("[FGG INT AVX] AVX Abort (PARAMS)\n");
      SE_int_split_SSE_dispatch(force, grid, q, spline, params, scale, bClearF);
      return;
    }
  // otherwise the preconditions for AVX codes are satisfied. 
  if(p==8){
    // specific for p=8
    __DISPATCHER_MSG("[FGG INT AVX] P=8\n");
    SE_int_split_AVX(force, grid, q, spline, params, scale, bClearF);
  }
}


#endif //GMX_X86_AVX_256
#endif // _SE_INT_AVX_256_SINGLE_
