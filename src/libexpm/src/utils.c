/* Copyright (C) 2013-2015 Drew Schmidt. All rights reserved.

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


/**
 * @file
 * @brief 
 * Copy matrix x ONTO y (y = x).
 *
 * @details
 * Simple wrapper around dlacpy.
 *
 * @param m
 * Number of rows of matrix x.
 * @param n
 * Number of cols of matrix x.
 * @param x
 * Input matrix.
 * @param y
 * Output matrix.
 */
void matcopy(int m, int n, double *x, double *y)
{
  char uplo = 'A';
  
  dlacpy_(&uplo, &m, &n, x, &m, y, &m);
}



/**
 * @file
 * @brief 
 * Construct identity matrix.
 *
 * @param n
 * Number of rows/cols of (square) matrix a.
 * @param a
 * Output matrix.
 */
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
// Products
// ---------------------------------------------------------

/**
 * @file
 * @brief 
 * Square matrix product.  Sets C = A*B
 * 
 * @details
 * Simple wrapper for dgemm for square matrices.
 *
 * @param n
 * Number of rows/cols of (square) matrices A, B, C.
 * @param A,B
 * Input matrices.
 * @param C
 * Output matrix.
 */
void matprod(int n, double *A, double *B, double *C)
{
  char trans = 'N';
  double one = 1.0, zero = 0.0;
  
  dgemm_(&trans, &trans, &n, &n, &n, &one, A, &n, B, &n, &zero, C, &n);
}



/**
 * @file
 * @brief 
 * y = A^pow * x
 * 
 * @details
 * Successively apply A (or A^T) to vector x, pow times.  dgemv is 
 * used for the products.
 *
 * @param trans
 * Input.  Determines if A or A^T is used.
 * @param pow
 * Power of matrix A.
 * @param n
 * Length of vectors x and y, and the number of rows/cols of square matrix A.
 * @param A
 * Input
 * @param x
 * 
 * @param y
 * Output
 */
void matvecprod(bool trans, int pow, int n, double *A, double *x, double *y)
{
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
    dgemv_(&transa, &n, &n, &one, A, &n, x, &ione, &zero, y, &ione);
    return;
  }
  
  
  tmp = malloc(n * sizeof(double));
  for (i=0; i<n; i++)
    y[i] = x[i];
  
  for (i=0; i<(pow-extra)/2; i+=2)
  {
    dgemv_(&transa, &n, &n, &one, A, &n, y, &ione, &zero, tmp, &ione);
    dgemv_(&transa, &n, &n, &one, A, &n, tmp, &ione, &zero, y, &ione);
  }
  
  
  if (extra)
    dgemv_(&transa, &n, &n, &one, A, &n, tmp, &ione, &zero, y, &ione);
  else
    matcopy(n, 1, tmp, y);
  
  
  free(tmp);
}










  // max(rowSums(abs(x)))
int vecnorm_inf(const int m, const int n, const double *x, double *norm, int *ind)
{
  int i, j;
  double *tmp;
  
  tmp = calloc(m, sizeof(*tmp));
  
  // tmp = rowSums(abs(x))
  for (j=0; j<n; j++)
  {
    for (i=0; i<m; i++)
      tmp[i] += fabs(x[i + j*m]);
  }
  
  *norm = tmp[0];
  *ind = 0;
  for (i=1; i<m; i++)
  {
    if (tmp[i] > *norm)
    {
      *norm = tmp[i];
      *ind = i;
    }
  }
  
  free(tmp);
  
  return 0;
}



