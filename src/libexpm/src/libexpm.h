#ifndef LIBEXPM_EXPM_H
#define LIBEXPM_EXPM_H


// expm.c
int libexpm_expm(const int p, const int n, double *x, double *ret);

// normest.c
double libexpm_normest_2_1(const int pow, const int n, double *A);
double libexpm_normest_2_4(const int pow, const int n, double *A);

// utils.c
int vecnorm_inf(const int m, const int n, const double *x, double *norm, int *ind);

#endif
