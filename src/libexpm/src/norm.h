#ifndef LIBEXPM_NORM_H
#define LIBEXPM_NORM_H


#define NORM_IDENTITY   1
#define NORM_ONE        2
#define NORM_FROBENIUS  3

#define NORM_ERROR_NORM   -1
#define NORM_ERROR_M      -2
#define NORM_ERROR_N      -3
#define NORM_ERROR_OOM    -2147483647

typedef int normtype_t;

double matnorm(normtype_t norm, const int m, const int n, const double *x, double *ret);


#endif
