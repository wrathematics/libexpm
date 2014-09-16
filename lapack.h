#ifndef EXPM_LAPACK_H
#define EXPM_LAPACK_H


#include "utils.h"


// BLAS
extern void F77_NAME(dgemm)(const char *transa, const char *transb, 
  const int *m, const int *n, const int *k, const double *alpha, 
  const double *a, const int *lda, const double *b, const int *ldb, 
  const double *beta, double *c, const int *ldc);

extern void F77_NAME(dscal)(const int *n, const double *alpha, double *dx, 
  const int *incx);


// LAPACK
extern void F77_NAME(dlacpy)(const char* uplo, const int* m, const int* n,
  const double* a, const int* lda, double* b, const int* ldb);

extern void F77_NAME(dgesv)(const int* n, const int* nrhs, double* a, 
  const int* lda, int* ipiv, double* b, const int* ldb, int* info);


#endif
