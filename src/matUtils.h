#ifndef MATUTILS_HEADER
#define MATUTILS_HEADER

#define CFree(x) if(x != NULL) free(x);

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Error.h>
#include <cmath>

void inverse(double** A, int N);
void **matrixMult(double **v1, int d11, int d12, double **v2, int d21, int d22, double **result);
double matrixDet(double **m, int dim);
void inverseR(double** mat, int N);
#endif
