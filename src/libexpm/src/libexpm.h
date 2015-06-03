#ifndef __LIBEXPM_H__
#define __LIBEXPM_H__


// restricted pointer double
typedef double* __restrict rpdouble_t;

// expm.c
int libexpm_expm_3_1(const int p, const int n, double *restrict x, double *restrict ret);
int libexpm_expm_5_1(const int p, const int n, double *restrict x, double *restrict ret);

// normest.c
double libexpm_normest_2_1(const int pow, const int n, double *A);
double libexpm_normest_2_4(const int pow, const int n, double *A);

// utils.c
int vecnorm_inf(const int m, const int n, const double *restrict x, double *restrict norm, int *restrict ind);

#endif
