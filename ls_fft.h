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

/*! \file ls_fft.h
 *  Interface for the LevelS FFT package.
 *
 *  \author Martin Reinecke
 */

#ifndef PLANCK_LS_FFT_H
#define PLANCK_LS_FFT_H

#include "c_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!\defgroup fftgroup FFT interface
This package is intended to calculate one-dimensional real or complex FFTs
with high accuracy and good efficiency even for lengths containing large
prime factors.

Before any FFT is executed, a plan must be generated for it. Plan creation
is designed to be fast, so that there is no significant overhead if the
plan is only used once or a few times.

The main component of the code is based on Paul N. Swarztrauber's FFTPACK in the
double precision incarnation by Hugh C. Pumphrey
(http://www.netlib.org/fftpack/dp.tgz).

I replaced the iterative sine and cosine calculations in radfg() and radbg()
by an exact calculation, which slightly improves the transform accuracy for
real FFTs with lengths containing large prime factors.

Since FFTPACK becomes quite slow for FFT lengths with large prime factors
(in the worst case of prime lengths it reaches \f$\mathcal{O}(n^2)\f$
complexity), I implemented Bluestein's algorithm, which computes a FFT of length
\f$n\f$ by several FFTs of length \f$n_2\ge 2n-1\f$ and a convolution. Since
\f$n_2\f$ can be chosen to be highly composite, this algorithm is more efficient
if \f$n\f$ has large prime factors. The longer FFTs themselves are then computed
using the FFTPACK routines.
Bluestein's algorithm was implemented according to the description on Wikipedia
(<a href="http://en.wikipedia.org/wiki/Bluestein%27s_FFT_algorithm">
http://en.wikipedia.org/wiki/Bluestein%27s_FFT_algorithm</a>).

\b Thread-safety:
All routines can be called concurrently; all information needed by
<tt>ls_fft</tt> is stored in the plan variable. However, using the same plan
variable on multiple threads simultaneously is not supported and will lead to
data corruption.
*/
/*! \{ */

typedef struct
  {
  double *work;
  size_t length, worksize;
  int bluestein;
  } complex_plan_i;

/*! The opaque handle type for complex-FFT plans. */
typedef complex_plan_i * complex_plan;

/*! Returns a plan for a complex FFT with \a length elements. */
complex_plan make_complex_plan (size_t length);
/*! Constructs a copy of \a plan. */
complex_plan copy_complex_plan (complex_plan plan);
/*! Destroys a plan for a complex FFT. */
void kill_complex_plan (complex_plan plan);
/*! Computes a complex forward FFT on \a data, using \a plan.
    \a Data has the form <tt>r0, i0, r1, i1, ...,
    r[length-1], i[length-1]</tt>. */
void complex_plan_forward (complex_plan plan, double *data);
/*! Computes a complex backward FFT on \a data, using \a plan.
    \a Data has the form <tt>r0, i0, r1, i1, ...,
    r[length-1], i[length-1]</tt>. */
void complex_plan_backward (complex_plan plan, double *data);

typedef struct
  {
  double *work;
  size_t length, worksize;
  int bluestein;
  } real_plan_i;

/*! The opaque handle type for real-FFT plans. */
typedef real_plan_i * real_plan;

/*! Returns a plan for a real FFT with \a length elements. */
real_plan make_real_plan (size_t length);
/*! Constructs a copy of \a plan. */
real_plan copy_real_plan (real_plan plan);
/*! Destroys a plan for a real FFT. */
void kill_real_plan (real_plan plan);
/*! Computes a real forward FFT on \a data, using \a plan
    and assuming the FFTPACK storage scheme:
    - on entry, \a data has the form <tt>r0, r1, ..., r[length-1]</tt>;
    - on exit, it has the form <tt>r0, r1, i1, r2, i2, ...</tt>
      (a total of \a length values). */
void real_plan_forward_fftpack (real_plan plan, double *data);
/*! Computes a real forward FFT on \a data, using \a plan
    and assuming the FFTPACK storage scheme:
    - on entry, \a data has the form <tt>r0, r1, i1, r2, i2, ...</tt>
    (a total of \a length values);
    - on exit, it has the form <tt>r0, r1, ..., r[length-1]</tt>. */
void real_plan_backward_fftpack (real_plan plan, double *data);
/*! Computes a real forward FFT on \a data, using \a plan
    and assuming the FFTW halfcomplex storage scheme:
    - on entry, \a data has the form <tt>r0, r1, ..., r[length-1]</tt>;
    - on exit, it has the form <tt>r0, r1, r2, ..., i2, i1</tt>. */
void real_plan_forward_fftw (real_plan plan, double *data);
/*! Computes a real backward FFT on \a data, using \a plan
    and assuming the FFTW halfcomplex storage scheme:
    - on entry, \a data has the form <tt>r0, r1, r2, ..., i2, i1</tt>.
    - on exit, it has the form <tt>r0, r1, ..., r[length-1]</tt>. */
void real_plan_backward_fftw (real_plan plan, double *data);
/*! Computes a real forward FFT on \a data, using \a plan
    and assuming a full-complex storage scheme:
    - on entry, \a data has the form <tt>r0, [ignored], r1, [ignored], ...,
      r[length-1], [ignored]</tt>;
    - on exit, it has the form <tt>r0, i0, r1, i1, ...,
      r[length-1], i[length-1]</tt>.
    */
void real_plan_forward_c (real_plan plan, double *data);
/*! Computes a real backward FFT on \a data, using \a plan
    and assuming a full-complex storage scheme:
    - on entry, \a data has the form <tt>r0, i0, r1, i1, ...,
      r[length-1], i[length-1]</tt>;
    - on exit, it has the form <tt>r0, 0, r1, 0, ..., r[length-1], 0</tt>. */
void real_plan_backward_c (real_plan plan, double *data);

void fftpack2halfcomplex (double *data, size_t n);
void halfcomplex2fftpack (double *data, size_t n);

/*! \} */

#ifdef __cplusplus
}
#endif

#endif
