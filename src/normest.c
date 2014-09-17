/* Copyright (C) 2014 Drew Schmidt. All rights reserved.

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


#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include "lapack.h"
#include "utils.h"

// Fortran codes for estimating the one-norm of a real or complex matrix, with applications to condition estimation (1988) 
// by Nicholas J. Higham 

// and

// A block algorithm for matrix 1-norm estimation with an application
// to 1-norm pseudospectra by Nicholas J. Higham and Francoise Tisseur


// Don't want to bring fortran into this (for the wrapper for ddot), and no real reason to use dgemv over dgemm
static inline double vecvecprod(int n, double *x, double *y)
{
  char transa = 'T', transb = 'F';
  int ione = 1;
  double one = 1.0, zero = 0.0;
  double dot;
  
  
  F77_NAME(dgemm)(&transa, &transb, &ione, &ione, &n, &one, x, &n, y, &n, &zero, &dot, &ione);
  
  return dot;
}



// Algorithm 2.1 from the Higham and Tisseur paper
#define ITMAX 9

double normest(const int pow, const int n, double *A)
{
  int i, k = 0;
  int ind;
  double tmp, norm;
  double ztx, znorm;
  double *x, *y, *s;
  
  
  if (pow == 1)
    return matnorm_1(A, n, n);
  
  x = malloc(n * sizeof(double));
  y = malloc(n * sizeof(double));
  s = malloc(n * sizeof(double));
  
  
  tmp = 1.0 / (double) n;
  for (i=0; i<n; i++)
    x[i] = tmp;
  
  while (k < ITMAX)
  {
    // y = A^pow x
    matvecprod(false, pow, n, A, x, y);
    
    norm = matnorm_1(y, n, 1);
    
    // s = sign(y)
    for (i=0; i<n; i++)
      s[i] = SGN(y[i]);
    
    // z = t(A)^pow s
    matvecprod(true, pow, n, A, s, y);
    
    znorm = vecnorm_inf(n, y, &ind)
    ztx = vecvecprod(n, y, x);
    
    if (znorm < ztx && k > 0)
      break;
    
    // x = e_j, where j is such that |z_j| = ||z||_inf
    for (i=0; i<n; i++)
      x[i] = 0.0;
    
    x[ind] = 1.0;
    
    k++;
  }
  
  free(x);
  free(y);
  free(s);
  
  return norm;
}





