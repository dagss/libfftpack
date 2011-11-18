/*
Copyright (c) 2004-2011, Max-Planck-Society
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fftpack.h"
#include "ls_fft.h"
#include "bluestein.h"

#define maxlen 4096

static void fill_random (double *data, double *odata, size_t length)
  {
  size_t m;
  for (m=0; m<length; ++m)
    data[m] = odata[m] = rand()/(RAND_MAX+1.0)-0.5;
  }
static void fill_random_c (double *data, double *odata, size_t length)
  {
  size_t m;
  for (m=0; m<length; ++m)
    {
    data[2*m] = odata[2*m] = rand()/(RAND_MAX+1.0)-0.5;
    data[2*m+1] = odata[2*m+1] = 0.;
    }
  }

static void normalize (double *data, size_t length, double norm)
  {
  size_t m;
  for (m=0; m<length; ++m)
    data[m] /= norm;
  }

static double errcalc (double *data, double *odata, size_t length)
  {
  size_t m;
  double sum, errsum;
  sum = errsum = 0;
  for (m=0; m<length; ++m)
    {
    errsum += (data[m]-odata[m])*(data[m]-odata[m]);
    sum += odata[m]*odata[m];
    }
  return sqrt(errsum/sum);
  }

static void test_real(void)
  {
  double data[2*maxlen], odata[2*maxlen];
  int length;
  double err;
  const double epsilon=3e-15;
  real_plan plan;
  for (length=1; length<=maxlen; ++length)
    {
    fill_random_c (data, odata, length);
    plan = make_real_plan (length);
    real_plan_forward_c(plan, data);
    real_plan_backward_c(plan, data);
    normalize (data, 2*length, length);
    kill_real_plan (plan);
    err = errcalc (data, odata, 2*length);
    if (err>epsilon) printf("problem at real length %i: %e\n",length,err);
    }
  for (length=1; length<=maxlen; ++length)
    {
    fill_random (data, odata, length);
    plan = make_real_plan (length);
    real_plan_forward_fftpack (plan, data);
    real_plan_backward_fftpack (plan, data);
    normalize (data, length, length);
    kill_real_plan (plan);
    err = errcalc (data, odata, length);
    if (err>epsilon) printf("problem at real length %i: %e\n",length,err);
    }
  for (length=1; length<=maxlen; ++length)
    {
    fill_random (data, odata, length);
    plan = make_real_plan (length);
    real_plan_forward_fftw (plan, data);
    real_plan_backward_fftw (plan, data);
    normalize (data, length, length);
    kill_real_plan (plan);
    err = errcalc (data, odata, length);
    if (err>epsilon) printf("problem at real length %i: %e\n",length,err);
    }
  for (length=1; length<=maxlen; ++length)
    {
    fill_random (data, odata, length);
    plan = make_real_plan (length);
    real_plan_forward_fftw (plan, data);
    halfcomplex2fftpack (data,plan->length);
    real_plan_backward_fftpack (plan, data);
    normalize (data, length, length);
    kill_real_plan (plan);
    err = errcalc (data, odata, length);
    if (err>epsilon) printf("problem at real length %i: %e\n",length,err);
    }
  for (length=1; length<=maxlen; ++length)
    {
    fill_random (data, odata, length);
    plan = make_real_plan (length);
    real_plan_forward_fftpack (plan, data);
    fftpack2halfcomplex (data,plan->length);
    real_plan_backward_fftw (plan, data);
    normalize (data, length, length);
    kill_real_plan (plan);
    err = errcalc (data, odata, length);
    if (err>epsilon) printf("problem at real length %i: %e\n",length,err);
    }
  }

static void test_complex(void)
  {
  double data[2*maxlen], odata[2*maxlen];
  int length;
  double err;
  const double epsilon=3e-15;
  complex_plan plan;
  for (length=1; length<=maxlen; ++length)
    {
    fill_random (data, odata, 2*length);
    plan = make_complex_plan (length);
    complex_plan_forward(plan, data);
    complex_plan_backward(plan, data);
    normalize (data, 2*length, length);
    kill_complex_plan (plan);
    err = errcalc (data, odata, 2*length);
    if (err>epsilon) printf("problem at complex length %i: %e\n",length,err);
    }
  }

int main(void)
  {
  test_real();
  test_complex();
  return 0;
  }
