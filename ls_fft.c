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

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "bluestein.h"
#include "fftpack.h"
#include "ls_fft.h"

complex_plan make_complex_plan (size_t length)
  {
  complex_plan plan = RALLOC(complex_plan_i,1);
  size_t pfsum = prime_factor_sum(length);
  double comp1 = (double)(length*pfsum);
  double comp2 = 2*3*length*log(3.*length);
  comp2*=3.; /* fudge factor that appears to give good overall performance */
  plan->length=length;
  plan->bluestein = (comp2<comp1);
  if (plan->bluestein)
    bluestein_i (length,&(plan->work),&(plan->worksize));
  else
    {
    plan->worksize=4*length+15;
    plan->work=RALLOC(double,4*length+15);
    cffti(length, plan->work);
    }
  return plan;
  }

complex_plan copy_complex_plan (complex_plan plan)
  {
  if (!plan) return NULL;
  {
  complex_plan newplan = RALLOC(complex_plan_i,1);
  *newplan = *plan;
  newplan->work=RALLOC(double,newplan->worksize);
  memcpy(newplan->work,plan->work,sizeof(double)*newplan->worksize);
  return newplan;
  }
  }

void kill_complex_plan (complex_plan plan)
  {
  DEALLOC(plan->work);
  DEALLOC(plan);
  }

void complex_plan_forward (complex_plan plan, double *data)
  {
  if (plan->bluestein)
    bluestein (plan->length, data, plan->work, -1);
  else
    cfftf (plan->length, data, plan->work);
  }

void complex_plan_backward (complex_plan plan, double *data)
  {
  if (plan->bluestein)
    bluestein (plan->length, data, plan->work, 1);
  else
    cfftb (plan->length, data, plan->work);
  }


real_plan make_real_plan (size_t length)
  {
  real_plan plan = RALLOC(real_plan_i,1);
  size_t pfsum = prime_factor_sum(length);
  double comp1 = .5*length*pfsum;
  double comp2 = 2*3*length*log(3.*length);
  comp2*=3; /* fudge factor that appears to give good overall performance */
  plan->length=length;
  plan->bluestein = (comp2<comp1);
  if (plan->bluestein)
    bluestein_i (length,&(plan->work),&(plan->worksize));
  else
    {
    plan->worksize=2*length+15;
    plan->work=RALLOC(double,2*length+15);
    rffti(length, plan->work);
    }
  return plan;
  }

real_plan copy_real_plan (real_plan plan)
  {
  if (!plan) return NULL;
  {
  real_plan newplan = RALLOC(real_plan_i,1);
  *newplan = *plan;
  newplan->work=RALLOC(double,newplan->worksize);
  memcpy(newplan->work,plan->work,sizeof(double)*newplan->worksize);
  return newplan;
  }
  }

void kill_real_plan (real_plan plan)
  {
  DEALLOC(plan->work);
  DEALLOC(plan);
  }

void real_plan_forward_fftpack (real_plan plan, double *data)
  {
  if (plan->bluestein)
    {
    size_t m;
    size_t n=plan->length;
    double *tmp = RALLOC(double,2*n);
    for (m=0; m<n; ++m)
      {
      tmp[2*m] = data[m];
      tmp[2*m+1] = 0.;
      }
    bluestein(n,tmp,plan->work,-1);
    data[0] = tmp[0];
    memcpy (data+1, tmp+2, (n-1)*sizeof(double));
    DEALLOC(tmp);
    }
  else
    rfftf (plan->length, data, plan->work);
  }

void fftpack2halfcomplex (double *data, size_t n)
  {
  size_t m;
  double *tmp = RALLOC(double,n);
  tmp[0]=data[0];
  for (m=1; m<(n+1)/2; ++m)
    {
    tmp[m]=data[2*m-1];
    tmp[n-m]=data[2*m];
    }
  if (!(n&1))
    tmp[n/2]=data[n-1];
  memcpy (data,tmp,n*sizeof(double));
  DEALLOC(tmp);
  }

void halfcomplex2fftpack (double *data, size_t n)
  {
  size_t m;
  double *tmp = RALLOC(double,n);
  tmp[0]=data[0];
  for (m=1; m<(n+1)/2; ++m)
    {
    tmp[2*m-1]=data[m];
    tmp[2*m]=data[n-m];
    }
  if (!(n&1))
    tmp[n-1]=data[n/2];
  memcpy (data,tmp,n*sizeof(double));
  DEALLOC(tmp);
  }

void real_plan_forward_fftw (real_plan plan, double *data)
  {
  real_plan_forward_fftpack (plan, data);
  fftpack2halfcomplex (data,plan->length);
  }

void real_plan_backward_fftpack (real_plan plan, double *data)
  {
  if (plan->bluestein)
    {
    size_t m;
    size_t n=plan->length;
    double *tmp = RALLOC(double,2*n);
    tmp[0]=data[0];
    tmp[1]=0.;
    memcpy (tmp+2,data+1, (n-1)*sizeof(double));
    if ((n&1)==0) tmp[n+1]=0.;
    for (m=2; m<n; m+=2)
      {
      tmp[2*n-m]=tmp[m];
      tmp[2*n-m+1]=-tmp[m+1];
      }
    bluestein (n, tmp, plan->work, 1);
    for (m=0; m<n; ++m)
      data[m] = tmp[2*m];
    DEALLOC(tmp);
    }
  else
    rfftb (plan->length, data, plan->work);
  }

void real_plan_backward_fftw (real_plan plan, double *data)
  {
  halfcomplex2fftpack (data,plan->length);
  real_plan_backward_fftpack (plan, data);
  }

void real_plan_forward_c (real_plan plan, double *data)
  {
  size_t m;
  size_t n=plan->length;

  if (plan->bluestein)
    {
    for (m=1; m<2*n; m+=2)
      data[m]=0;
    bluestein (plan->length, data, plan->work, -1);
    data[1]=0;
    for (m=2; m<n; m+=2)
      {
      double avg;
      avg = 0.5*(data[2*n-m]+data[m]);
      data[2*n-m] = data[m] = avg;
      avg = 0.5*(data[2*n-m+1]-data[m+1]);
      data[2*n-m+1] = avg;
      data[m+1] = -avg;
      }
    if ((n&1)==0) data[n+1] = 0.;
    }
  else
    {
/* using "m+m" instead of "2*m" to avoid a nasty bug in Intel's compiler */
    for (m=0; m<n; ++m) data[m+1] = data[m+m];
    rfftf (n, data+1, plan->work);
    data[0] = data[1];
    data[1] = 0;
    for (m=2; m<n; m+=2)
      {
      data[2*n-m]   =  data[m];
      data[2*n-m+1] = -data[m+1];
      }
    if ((n&1)==0) data[n+1] = 0.;
    }
  }

void real_plan_backward_c (real_plan plan, double *data)
  {
  size_t n=plan->length;

  if (plan->bluestein)
    {
    size_t m;
    data[1]=0;
    for (m=2; m<n; m+=2)
      {
      double avg;
      avg = 0.5*(data[2*n-m]+data[m]);
      data[2*n-m] = data[m] = avg;
      avg = 0.5*(data[2*n-m+1]-data[m+1]);
      data[2*n-m+1] = avg;
      data[m+1] = -avg;
      }
    if ((n&1)==0) data[n+1] = 0.;
    bluestein (plan->length, data, plan->work, 1);
    for (m=1; m<2*n; m+=2)
      data[m]=0;
    }
  else
    {
    ptrdiff_t m;
    data[1] = data[0];
    rfftb (n, data+1, plan->work);
    for (m=n-1; m>=0; --m)
      {
      data[2*m]   = data[m+1];
      data[2*m+1] = 0.;
      }
    }
  }
