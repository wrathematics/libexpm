#ifndef LIBEXPM_UTILS_H
#define LIBEXPM_UTILS_H


#include <stdbool.h>


#define SGN(x) (x>=0?1:-1)
#define SGNEXP(x,pow) (x==0?(pow==0?1:0):(x>0?1:(pow%2==0?1:(-1))))

#define MIN(a,b) (a<b?a:b)

void matcopy(int n, double *A, double *B);
void mateye(const unsigned int n, double *a);
double matnorm_1(const double *x, const int m, const int n);
void matprod(int n, double *a, double *b, double *c);
void matvecprod(bool trans, int pow, int n, double *a, double *b, double *c);


#endif
