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


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "libexpm.h"
#include "utils.h"
#include "lapack.h"
#include "norm.h"


#define NTHETA 5

static int matexp_scale_factor(const double *x, const int n)
{
  int i;
  double x_1;
  const double theta[] = {1.5e-2, 2.5e-1, 9.5e-1, 2.1e0, 5.4e0};
  
  matnorm(NORM_ONE, n, n, x, &x_1);
  
  for (i=0; i<NTHETA; i++)
  {
    if (x_1 <= theta[i])
      return 0;
  }
  
  i = (int) ceil(log2(x_1/theta[4]));
  
  return 1 << i;
}



// Matrix power by squaring: P = A^b (A is garbage on exit)
static void matpow_by_squaring(double *A, int n, int b, double *P)
{
  double *TMP;
  
  mateye(n, P);
  
  // Trivial cases
  if (b == 0)
    return;
  
  if (b == 1)
  {
    matcopy(n, n, A, P);
    return;
  }
  
  
  // General case
  TMP = malloc(n*n*sizeof(double));
  
  while (b)
  {
    if (b&1)
    {
      matprod(n, P, A, TMP);
      matcopy(n, n, TMP, P);
    }
    
    b >>= 1;
    matprod(n, A, A, TMP);
    matcopy(n, n, TMP, A);
  }
  
  free(TMP);
}



// -------------------------------------------------------- 
// Matrix Exponentiation via Pade' Approximations
// -------------------------------------------------------- 

const double matexp_pade_coefs[14] =
{
  1.0,
  0.5,
  0.12,
  1.833333333333333333333e-2,
  1.992753623188405797101e-3,
  1.630434782608695652174e-4,
  1.035196687370600414079e-5,
  5.175983436853002070393e-7,
  2.043151356652500817261e-8,
  6.306022705717595115002e-10,
  1.483770048404140027059e-11,
  2.529153491597965955215e-13,
  2.810170546219962172461e-15,
  1.544049750670308885967e-17
};



/* r_m(x) = p_m(x) / q_m(x), where
   p_m(x) = sum_{j=0}^m (2m-j)!m!/(2m)!/(m-j)!/j! * x^j
  
   and q_m(x) = p_m(-x)
*/


// Exponentiation via Pade' expansion
static void matexp_pade(const int n, const int p, double *restrict A, double *restrict N)
{
  int i, j, info = 0;
  int *ipiv;
  int sgn;
  double tmp, tmpj;
  double *B, *C, *D;
  
  if (n <= 0) return;
  
  // Power of A
  B = calloc(n*n, sizeof(*B));
  assert(B != NULL);
  
  // Temporary storage for matrix multiplication
  C = malloc(n*n * sizeof(*C));
  assert(C != NULL);
  
  D = calloc(n*n, sizeof(*D));
  assert(D != NULL);
  
  matcopy(n, n, A, C);
  
  memset(D, 0, n*n*sizeof(*D));
  
  
  i = 0;
  while (i < n*n)
  {
    N[i] = 1.0;
    D[i] = 1.0;
    
    i += n+1;
  }
  
  
  // Fill N and D
  for (i=1; i<=p; i++)
  {
    // C = A*B
    if (i > 1)
      matprod(n, A, B, C);
    
    // Update matrices
    sgn = SGNEXP(-1, i);
    tmp = matexp_pade_coefs[i];
    
    /* Performs the following actions:
        B = C
        N = pade_coef[i] * C
        D = (-1)^j * pade_coef[i] * C
    */
    for (j=0; j<n*n; j++)
    {
      tmpj = C[j];
      B[j] = tmpj;
      
      tmpj *= tmp;
      
      N[j] += tmpj;
      D[j] += sgn*tmpj;
    }
  }
  
  // R <- inverse(D) %*% N
  ipiv = calloc(n, sizeof(*ipiv));
  assert(ipiv != NULL);
  
  dgesv_(&n, &n, D, &n, ipiv, N, &n, &info);
  
  
  free(B);
  free(C);
  free(D);
  free(ipiv);
}



/**
 * @file
 * @brief 
 * Matrix exponentiation via scaling and squaring.
 * 
 * @description
 * Algorithm 3.1 from the paper "A New Scaling and Squaring Algorithm
 * for the Matrix Exponential" by Awad H. Al-Mohy and Nicholas J. Higham.
 * 
 * @param n
 * Number of rows/cols of (square) matrix x.
 * @param p
 * Order of the Pade' approximation.  Must satisfy 0<p<=13.  
 * Useful values are p=6, 8, 12.
 * @param t
 * Scaling factor for x (t=1 canonical).
 * @param x
 * Modified Input (square) matrix, stored in column-major format.
 * On successful function return, the values of x are garbage.
 * @param ret
 * On successful return, ret = expm(x).
 *
 * @return
 * The return value indicates the status of the function.
 */
int libexpm_expm_3_1(const int p, const int n, double *restrict x, double *restrict ret)
{
  int m;
  int nn = n*n;
  int one = 1;
  double tmp;
  
  m = matexp_scale_factor(x, n);
  
  if (m == 0)
  {
    matexp_pade(n, p, x, ret);
    return 0;
  }
  
  tmp = 1. / ((double) m);
  dscal_(&nn, &tmp, x, &one);
  
  
  matexp_pade(n, p, x, ret);
  
  matcopy(n, n, ret, x);
  
  matpow_by_squaring(x, n, m, ret);
  
  return 0;
}





/**
 * @file
 * @brief 
 * Matrix exponentiation via new scaling and squaring algorithm.
 * 
 * @description
 * Algorithm 5.1 from the paper "A New Scaling and Squaring Algorithm
 * for the Matrix Exponential" by Awad H. Al-Mohy and Nicholas J. Higham.
 * 
 * @param n
 * Number of rows/cols of (square) matrix x.
 * @param p
 * Order of the Pade' approximation.  Must satisfy 0<p<=13.  
 * Useful values are p=6, 8, 12.
 * @param t
 * Scaling factor for x (t=1 canonical).
 * @param x
 * Modified Input (square) matrix, stored in column-major format.
 * On successful function return, the values of x are garbage.
 * @param ret
 * On successful return, ret = expm(x).
 *
 * @return
 * The return value indicates the status of the function.
 */
int libexpm_expm_5_1(const int p, const int n, double *restrict x, double *restrict ret)
{
  int m;
  int nn = n*n;
  int one = 1;
  double tmp;
  
  m = matexp_scale_factor(x, n);
  
  if (m == 0)
  {
    matexp_pade(n, p, x, ret);
    return 0;
  }
  
  tmp = 1. / ((double) m);
  dscal_(&nn, &tmp, x, &one);
  
  
  matexp_pade(n, p, x, ret);
  
  matcopy(n, n, ret, x);
  
  matpow_by_squaring(x, n, m, ret);
  
  return 0;
}

