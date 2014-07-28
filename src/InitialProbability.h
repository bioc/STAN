#ifndef INITPROB_HEADER
#define INITPROB_HEADER

#include <stdio.h>
#include <stdlib.h>
#include "DebugConstants.h"
#include <R.h>
#include <Rdefines.h>
#include "RAccessUtils.h"

using namespace std;

class InitialProbability
{
protected:
    int K;
    double* updateNumeratorPI;
    double* pi;
    
public:

    InitialProbability(double *pi, int K);
    virtual ~InitialProbability();

    double* getInitialProb();
    void updateSample(double** gamma, int i);
    void updateSampleCoupled(double** gamma, int i, int* couples, SEXP bidirOptimParams);
    void update(int nsample, SEXP bidirOptimParams);
    int getK();
    void finalize();
};
#endif
