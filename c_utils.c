/*
Copyright (c) 2008-2011, Max-Planck-Society
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
 *  libc_utils is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Convenience functions
 *
 *  Author: Martin Reinecke
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "c_utils.h"

void util_fail_ (const char *file, int line, const char *func, const char *msg)
  {
  fprintf(stderr,"%s, %i (%s):\n%s\n",file,line,func,msg);
  exit(1);
  }
void util_warn_ (const char *file, int line, const char *func, const char *msg)
  {
  fprintf(stderr,"%s, %i (%s):\n%s\n",file,line,func,msg);
  exit(1);
  }

void *util_malloc_ (size_t sz)
  {
  void *res;
  if (sz==0) return NULL;
  res = malloc(sz);
  UTIL_ASSERT(res,"malloc() failed");
  return res;
  }
void util_free_ (void *ptr)
  { if ((ptr)!=NULL) free(ptr); }
