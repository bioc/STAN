
#ifndef MEMORYALLOCATION_HEADER
#define MEMORYALLOCATION_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <R_ext/Error.h>
#include "ParamContainerEmissions.h"

double** allocateNumericMatrix(int d1, int d2);
double* allocateNumericVector(int d);

int allocateMemAlpha(double*** alpha, int maxLen, int K);
int allocateMemBeta(double*** beta, int maxLen, int K);
int allocateMemRescFac(double** c, int maxLen, int K);
int allocateMemGamma(double*** gamma, int maxLen, int K);
int allocateMemXsi(double**** xsi, int maxLen, int K);
int allocateMemEmissionProb(double*** emissionProb, int maxLen, int K);

int deallocateMemAlpha(double** alpha, int maxLen, int K);
int deallocateMemBeta(double** beta, int maxLen, int K);
int deallocateMemRescFac(double* c, int maxLen, int K);
int deallocateMemGamma(double** gamma, int maxLen, int K);
int deallocateMemXsi(double*** xsi, int maxLen, int K);
int deallocateMemEmissionProb(double** emissionProb, int maxLen, int K);
#endif
