/*
Copyright (c) 2010-2011, Max-Planck-Society
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
  fftpack.h : function declarations for fftpack.c
  Algorithmically based on Fortran-77 FFTPACK by Paul N. Swarztrauber
  (Version 4, 1985).

  C port by Martin Reinecke (2010)
 */

#ifndef PLANCK_FFTPACK_H
#define PLANCK_FFTPACK_H

#include "c_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! forward complex transform */
void cfftf(size_t N, double complex_data[], double wrk[]);
/*! backward complex transform */
void cfftb(size_t N, double complex_data[], double wrk[]);
/*! initializer for complex transforms */
void cffti(size_t N, double wrk[]);

/*! forward real transform */
void rfftf(size_t N, double data[], double wrk[]);
/*! backward real transform */
void rfftb(size_t N, double data[], double wrk[]);
/*! initializer for real transforms */
void rffti(size_t N, double wrk[]);

#ifdef __cplusplus
}
#endif

#endif
