#ifndef _SE_INT_AVX_256_SINGLE_
#define _SE_INT_AVX_256_SINGLE_

#if GMX_DOUBLE==0
/* SE AVX 256 single integration */
#if GMX_SIMD_X86_AVX_256
#include "se.h"

#include "gromacs/simd/simd.h"
static inline float gmx_simdcall
reduce(__m256 a)
{
    __m128 t0;
    t0 = _mm_add_ps(_mm256_castps256_ps128(a), _mm256_extractf128_ps(a, 0x1));
    t0 = _mm_add_ps(t0, _mm_permute_ps(t0, _MM_SHUFFLE(1, 0, 3, 2)));
    t0 = _mm_add_ss(t0, _mm_permute_ps(t0, _MM_SHUFFLE(0, 3, 2, 1)));
    return *reinterpret_cast<float *>(&t0);
}

// -----------------------------------------------------------------------------
static void SE_int_split_AVX(rvec * gmx_restrict force,  real * gmx_restrict grid, 
			     real * gmx_restrict q,
			     splinedata_t * gmx_restrict spline,
			     const SE_FGG_params* gmx_restrict params, 
			     real scale, gmx_bool bClearF,
			     const pme_atomcomm_t * gmx_restrict atc,
			     const gmx_pme_t * gmx_restrict pme)
{
  // unpack params
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
  int * idxptr;
  int i0,j0,k0;
  int index_x, index_xy;

  float qm, h3=h*h*h;

  int pny   = pme->pmegrid_ny;
  int pnz   = pme->pmegrid_nz;

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

  for(m=0; m<N; m++){
    mm  = spline->ind[m];
    qm = q[mm];
    idxptr = atc->idx[mm];
    i0 = idxptr[XX];
    j0 = idxptr[YY];
    k0 = idxptr[ZZ];
    idx = i0*pny*pnz+j0*pnz+k0;

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
    rZFZ0 = _mm256_load_ps(zfz + m*8   );

    if(idx%8==0){ // grid[idx] is 32-aligned so vectorization simple
      for(i = 0; i<8; i++){
	index_x = (i0+i)*pny*pnz;
	for(j = 0; j<8; j++){
	  index_xy = index_x + (j0+j)*pnz;
	  rC  = _mm256_set1_ps( zx[m*8+i]*zy[m*8 + j]);
	  rCX = _mm256_set1_ps( zx[m*8+i]*zy[m*8 + j] * zfx[m*8 + i]);
	  rCY = _mm256_set1_ps( zx[m*8+i]*zy[m*8 + j] * zfy[m*8 + j]);

	  rH0  = _mm256_load_ps( grid+index_xy +k0  );
	  rZS0 = _mm256_load_ps( zs + idx_zs     );

	  __m256 fac = _mm256_mul_ps(rH0, _mm256_mul_ps(rZZ0, rZS0));

	  rFX =_mm256_add_ps(rFX,_mm256_mul_ps(fac,rCX)); 
	  rFY =_mm256_add_ps(rFY,_mm256_mul_ps(fac,rCY));
	  rFZ =_mm256_add_ps(rFZ,_mm256_mul_ps(fac,_mm256_mul_ps(rC,rZFZ0))); 

#ifdef CALC_ENERGY
	  rP =_mm256_add_ps(rP,_mm256_mul_ps(fac,rC));
#endif
	  idx_zs +=8;
	}
      }
    }
    else{ // grid[idx] not 32-aligned, so use non-aligned loads
      for(i = 0; i<8; i++){
	index_x = (i0+i)*pny*pnz;
	for(j = 0; j<8; j++){
	  index_xy = index_x + (j0+j)*pnz;
	  rC  = _mm256_set1_ps( zx[m*8+i]*zy[m*8 + j]);
	  rCX = _mm256_set1_ps( zx[m*8+i]*zy[m*8 + j] * zfx[m*8 + i]);
	  rCY = _mm256_set1_ps( zx[m*8+i]*zy[m*8 + j] * zfy[m*8 + j]);
#ifdef CALC_ENERGY
	  rCP = _mm256_set1_ps( zx[m*8+i]*zy[m*8 + j]);
#endif

	  rH0  = _mm256_loadu_ps( grid+index_xy + k0 );
	  rZS0 = _mm256_load_ps( zs + idx_zs);
		 		    
	  __m256 fac = _mm256_mul_ps(rH0, _mm256_mul_ps(rZZ0, rZS0));

	  rFX =_mm256_add_ps(rFX,_mm256_mul_ps(fac,rCX)); 
	  rFY =_mm256_add_ps(rFY,_mm256_mul_ps(fac,rCY));
	  rFZ =_mm256_add_ps(rFZ,_mm256_mul_ps(fac,_mm256_mul_ps(rC,rZFZ0))); 
	 
#ifdef CALC_ENERGY
	  rP =_mm256_add_ps(rP,_mm256_mul_ps(fac,rC));
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

    force[m][XX] += -qm*scale*h3*(reduce(rFX));
    force[m][YY] += -qm*scale*h3*(reduce(rFY));
    force[m][ZZ] += -qm*scale*h3*(reduce(rFZ));

#ifdef CALC_ENERGY
    _mm256_store_ps(s,rP);
    st->phi[m] = -scale*h3*(reduce(rP));
#endif	
  }
}

// ----------------------------------------------------------------------------
static void 
SE_int_split_AVX_dispatch(rvec* force, real* grid, real* q,
			  splinedata_t *spline,
			  const SE_FGG_params* params, real scale,
			  gmx_bool bClearF,
			  const pme_atomcomm_t *atc,
			  const gmx_pme_t *pme)
{
  const int p = params->P;
  const int incrj = params->dims[2]; // middle increment
  const int incri = params->npdims[2]*(params->dims[1]);// outer increment

  // if P, incri or increments are not divisible by 4, fall back on vanilla
  if( isnot_div_by_8(p) || isnot_div_by_8(incri) || isnot_div_by_8(incrj) || (p%8)!=0)
    {
      __DISPATCHER_MSG("[FGG INT AVX SINGLE] AVX Abort (PARAMS)\n");
      SE_int_split_SSE_dispatch(force, grid, q, spline, params, scale, 
				bClearF, atc, pme);
      return;
    }
  // otherwise the preconditions for AVX codes are satisfied. 
  if(p==8){
    // specific for p=8
    __DISPATCHER_MSG("[FGG INT AVX SINGLE] P=8\n");
    //SE_int_split_SSE_dispatch(force, grid, q, spline, params, scale, bClearF);
    SE_int_split_AVX(force, grid, q, spline, params, scale, bClearF,
		     atc, pme);
  }
}

#endif // not GMX_DOUBLE
#endif // GMX_SIMD_X86_AVX_256
#endif // _SE_INT_AVX_256_SINGLE_
