#ifndef EXPM_LAPACK_H
#define EXPM_LAPACK_H


void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
void dlacpy_(char *uplo, int *m, int *n, double *a, int *lda, double *b, int *ldb);
void dscal_(int *n, double *a, double *x, int *incx);
void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);


#endif
