/*
Copyright (c) 2005-2011, Max-Planck-Society
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the Max-Planck-Society nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE MAX-PLANCK-SOCIETY BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
 *  libfftpack is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  \author Martin Reinecke
 */

#include <math.h>
#include <stdlib.h>
#include "fftpack.h"
#include "bluestein.h"

/* returns the sum of all prime factors of n */
size_t prime_factor_sum (size_t n)
  {
  size_t result=0,x,limit,tmp;
  while (((tmp=(n>>1))<<1)==n)
    { result+=2; n=tmp; }

  limit=(size_t)sqrt(n+0.01);
  for (x=3; x<=limit; x+=2)
  while ((tmp=(n/x))*x==n)
    {
    result+=x;
    n=tmp;
    limit=(size_t)sqrt(n+0.01);
    }
  if (n>1) result+=n;

  return result;
  }

/* returns the smallest composite of 2, 3 and 5 which is >= n */
static size_t good_size(size_t n)
  {
  size_t f2, f23, f235, bestfac=2*n;
  if (n<=6) return n;

  for (f2=1; f2<bestfac; f2*=2)
    for (f23=f2; f23<bestfac; f23*=3)
      for (f235=f23; f235<bestfac; f235*=5)
        if (f235>=n) bestfac=f235;
  return bestfac;
  }

void bluestein_i (size_t n, double **tstorage, size_t *worksize)
  {
  static const double pi=3.14159265358979323846;
  size_t n2=good_size(n*2-1);
  size_t m, coeff;
  double angle, xn2;
  double *bk, *bkf, *work;
  double pibyn=pi/n;
  *worksize=2+2*n+8*n2+16;
  *tstorage = RALLOC(double,2+2*n+8*n2+16);
  ((size_t *)(*tstorage))[0]=n2;
  bk  = *tstorage+2;
  bkf = *tstorage+2+2*n;
  work= *tstorage+2+2*(n+n2);

/* initialize b_k */
  bk[0] = 1;
  bk[1] = 0;

  coeff=0;
  for (m=1; m<n; ++m)
    {
    coeff+=2*m-1;
    if (coeff>=2*n) coeff-=2*n;
    angle = pibyn*coeff;
    bk[2*m] = cos(angle);
    bk[2*m+1] = sin(angle);
    }

/* initialize the zero-padded, Fourier transformed b_k. Add normalisation. */
  xn2 = 1./n2;
  bkf[0] = bk[0]*xn2;
  bkf[1] = bk[1]*xn2;
  for (m=2; m<2*n; m+=2)
    {
    bkf[m]   = bkf[2*n2-m]   = bk[m]   *xn2;
    bkf[m+1] = bkf[2*n2-m+1] = bk[m+1] *xn2;
    }
  for (m=2*n;m<=(2*n2-2*n+1);++m)
    bkf[m]=0.;
  cffti (n2,work);
  cfftf (n2,bkf,work);
  }

void bluestein (size_t n, double *data, double *tstorage, int isign)
  {
  size_t n2=*((size_t *)tstorage);
  size_t m;
  double *bk, *bkf, *akf, *work;
  bk  = tstorage+2;
  bkf = tstorage+2+2*n;
  work= tstorage+2+2*(n+n2);
  akf = tstorage+2+2*n+6*n2+16;

/* initialize a_k and FFT it */
  if (isign>0)
    for (m=0; m<2*n; m+=2)
      {
      akf[m]   = data[m]*bk[m]   - data[m+1]*bk[m+1];
      akf[m+1] = data[m]*bk[m+1] + data[m+1]*bk[m];
      }
  else
    for (m=0; m<2*n; m+=2)
      {
      akf[m]   = data[m]*bk[m]   + data[m+1]*bk[m+1];
      akf[m+1] =-data[m]*bk[m+1] + data[m+1]*bk[m];
      }
  for (m=2*n; m<2*n2; ++m)
    akf[m]=0;

  cfftf (n2,akf,work);

/* do the convolution */
  if (isign>0)
    for (m=0; m<2*n2; m+=2)
      {
      double im = -akf[m]*bkf[m+1] + akf[m+1]*bkf[m];
      akf[m  ]  =  akf[m]*bkf[m]   + akf[m+1]*bkf[m+1];
      akf[m+1]  = im;
      }
  else
    for (m=0; m<2*n2; m+=2)
      {
      double im = akf[m]*bkf[m+1] + akf[m+1]*bkf[m];
      akf[m  ]  = akf[m]*bkf[m]   - akf[m+1]*bkf[m+1];
      akf[m+1]  = im;
      }


/* inverse FFT */
  cfftb (n2,akf,work);

/* multiply by b_k* */
  if (isign>0)
    for (m=0; m<2*n; m+=2)
      {
      data[m]   = bk[m]  *akf[m] - bk[m+1]*akf[m+1];
      data[m+1] = bk[m+1]*akf[m] + bk[m]  *akf[m+1];
      }
  else
    for (m=0; m<2*n; m+=2)
      {
      data[m]   = bk[m]  *akf[m] + bk[m+1]*akf[m+1];
      data[m+1] =-bk[m+1]*akf[m] + bk[m]  *akf[m+1];
      }
  }
