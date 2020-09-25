#ifndef SE_KAISER_FFT_H
#define SE_KAISER_FFT_H

#include "gromacs/fft/fft.h"
#include "fstream"

// Do 1d FFT
void compute_fft1d(real *fftgrid, t_complex *cfftgrid, int M)
{
  gmx_fft_t pfft;
  gmx_fft_init_1d_real(&pfft, M, GMX_FFT_FLAG_NONE);

  gmx_fft_1d_real(pfft, GMX_FFT_REAL_TO_COMPLEX,
		  (void*) fftgrid, (void*) cfftgrid);

}

void compute_fft_kaiser(real *cfftgrid, SE_opt opt)
{
  real L = opt.box[XX];
  int  M = opt.M;  // FIXME: M is equal in all directions
  real h = L/M;
  int  P = opt.P;
  real w = P/2.;
  const real ow2   = 1./(w*w);
  real beta = opt.beta;
  real *fft1d;
  t_complex  *cfft1d;

  snew(fft1d, M);
  snew(cfft1d, M);

  // Construct 1D real kaiser window in a 3D vector
  for(int ix=0; ix<M; ix++){
    real x = 0+h*ix;
    real v = (x-L/2.)/h;
    fft1d[ix] = kaiser(v, ow2, beta)*h;
  }
  
  /* FIXME: Assuming that M is even, the imaginary part of
   * cfft1d is zero or below machine precision. Therefore,
   * we use a real to real fft.
   */
  compute_fft1d(fft1d, cfft1d, M);
  for(int ix=0; ix<M/2; ix++)
    cfft1d[M-ix].re = cfft1d[ix].re;

  for(int ix=0; ix<M; ix++)
    cfftgrid[ix] = cfft1d[ix].re;
  /* for(int ix=0; ix<M; ix++) */
  /*   printf("KAISER %.6g\t %.6g\n",fft1d[ix],cfft1d[ix].re); */

  /* for(int iy=0; iy<M; iy++) */
  /*   for(int iz=0; iz<M; iz++) */
  /*     for(int ix=0; ix<M; ix++) */
  /* 	{ */
  /* 	  int fftidx = ix*M*M+iy*M+iz; */
  /* 	  real v = */
  /* 	    cfft1d[ix].re*cfft1d[ix].re* */
  /* 	    cfft1d[iy].re*cfft1d[iy].re* */
  /* 	    cfft1d[iz].re*cfft1d[iz].re; */
	  
  /* 	  //  	if(iz!=(M/2-1)) */
  /* 	  cfftgrid[fftidx] = 1/v; */
  /* 	} */

  /* std::ofstream fp; */
  /* fp.open("KAISER.txt", std::ios_base::app); */
  
  /* for (int ix=0; ix<M*M*M; ix++) */
  /*   fp << cfftgrid[ix] << "\n"; */

  /* fp.close(); */
  
  //sfree(fft1d);
  //sfree(cfft1d);
}

#endif // SE_KAISER_FFT_H
