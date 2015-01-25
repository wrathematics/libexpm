#ifndef LIBEXPM_LAPACK_H
#define LIBEXPM_LAPACK_H


// FIXME
#define F77_NAME(x) x ## _

// BLAS
void F77_NAME(dgemm)(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
void F77_NAME(dscal)(int *n, double *a, double *x, int *incx);

// LAPACK
void F77_NAME(dlacpy)(char *uplo, int *m, int *n, double *a, int *lda, double *b, int *ldb);
void F77_NAME(dgesv)(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);


#endif
