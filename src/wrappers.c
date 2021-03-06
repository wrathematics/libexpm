/* Copyright (C) 2014 Drew Schmidt. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */


#include <R.h>
#include <Rinternals.h>
#include "libexpm/src/libexpm.h"


SEXP R_libexpm_normest(SEXP x, SEXP pow)
{
  const int n = nrows(x);
  int i;
  SEXP ret;
  PROTECT(ret = allocVector(REALSXP, 1));
  
  REAL(ret)[0] = libexpm_normest_2_4(INTEGER(pow)[0], n, REAL(x));
  
  UNPROTECT(1);
  return ret;
}



SEXP R_libexpm_expm_3_1(SEXP x, SEXP p)
{
  const int n = nrows(x);
  int i;
  double *x_cp;
  SEXP R;
  
  PROTECT(R = allocMatrix(REALSXP, n, n));
  
  x_cp = malloc(n*n * sizeof(*x_cp));
  memcpy(x_cp, REAL(x), n*n*sizeof(*x_cp));
  
  libexpm_expm_3_1(INTEGER(p)[0], n, x_cp, REAL(R));
  
  free(x_cp);
  
  UNPROTECT(1);
  return R;
}




SEXP R_libexpm_expm_5_1(SEXP x, SEXP p)
{
  const int n = nrows(x);
  int i;
  double *x_cp;
  SEXP R;
  
  PROTECT(R = allocMatrix(REALSXP, n, n));
  
  x_cp = malloc(n*n * sizeof(*x_cp));
  memcpy(x_cp, REAL(x), n*n*sizeof(*x_cp));
  
  libexpm_expm_5_1(INTEGER(p)[0], n, x_cp, REAL(R));
  
  free(x_cp);
  
  UNPROTECT(1);
  return R;
}

