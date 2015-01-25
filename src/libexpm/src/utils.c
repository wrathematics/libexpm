/* Copyright (C) 2013-2014 Drew Schmidt. All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
  
    * Redistributions of source code must retain the above copyright notice, 
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */


/* Matrix exponentiation algorithm from:
   "New Scaling and Squaring Algorithm for the Matrix Exponential", by
   Awad H. Al-Mohy and Nicholas J. Higham, August 2009
*/


#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "lapack.h"
#include "utils.h"


// ---------------------------------------------------------
// Misc
// ---------------------------------------------------------

// Copy A ONTO B, i.e. B = A
void matcopy(int n, double *A, double *B)
{
  char uplo = 'A';
  
  
  F77_NAME(dlacpy)(&uplo, &n, &n, A, &n, B, &n);
}



// Identity matrix
void mateye(const unsigned int n, double *a)
{
  int i;
  
  
  for (i=0; i<n*n; i++)
    a[i] = 0.0;
  
  i = 0;
  while (i < n*n)
  {
    a[i] = 1.0;
    i += n+1;
  }
}



// ---------------------------------------------------------
// Explicit 1-Norms
// ---------------------------------------------------------


// Full 1-norm
double matnorm_1(const double *x, const int m, const int n)
{
  int i, j;
  double norm = -1.;
  double tmp;
  
  
  // max(colSums(abs(x))) 
  for (j=0; j<n; j++)
  {
    tmp = 0;
    
    for (i=0; i<m; i++)
      tmp += fabs(x[i + j*m]);
    
    if (tmp > norm)
      norm = tmp;
  }
  
  return norm;
}



double vecnorm_inf(const int n, const double *x, int *ind)
{
  int i;
  double tmp, norm = -1.;
  
  
  // max(abs(x))
  for (i=0; i<n; i++)
  {
    tmp = fabs(x[i]);
    
    if (tmp > norm)
    {
      norm = tmp;
      *ind = i;
    }
  }
  
  return norm;
}



// ---------------------------------------------------------
// Products
// ---------------------------------------------------------

// C = A * B for square matrices
void matprod(int n, double *a, double *b, double *c)
{
  char trans = 'N';
  double one = 1.0, zero = 0.0;
  
  
  F77_NAME(dgemm)(&trans, &trans, &n, &n, &n, &one, a, &n, b, &n, &zero, c, &n);
}



// Matrix a times vector b
void matvecprod(bool trans, int pow, int n, double *a, double *b, double *c)
{
  char notrans = 'N';
  char transa;
  int i;
  int ione = 1;
  int extra = pow % 2;
  double one = 1.0, zero = 0.0;
  double *tmp;
  
  
  if (transa)
    transa = 'T';
  else
    transa = 'N';
  
  if (pow == 1)
  {
    F77_NAME(dgemm)(&transa, &notrans, &n, &ione, &n, &one, a, &n, b, &n, &zero, c, &n);
    return;
  }
  
  
  tmp = malloc(n*sizeof(*tmp));
  
  for (i=0; i<pow-extra; i++)
  {
    F77_NAME(dgemm)(&transa, &notrans, &n, &ione, &n, &one, a, &n, tmp, &n, &zero, c, &n);
    F77_NAME(dgemm)(&transa, &notrans, &n, &ione, &n, &one, a, &n, c, &n, &zero, tmp, &n);
  }
  
  if (extra)
    F77_NAME(dgemm)(&transa, &notrans, &n, &ione, &n, &one, a, &n, tmp, &n, &zero, c, &n);
  else
    matcopy(n, tmp, c);
  
  free(tmp);
}


