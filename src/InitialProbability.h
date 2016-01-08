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
        double* pi;                               //!<@brief K-dimensional vector containing initial state probabilities.
// TODO: also constraints here!

    public:

        InitialProbability(double *pi, int K);

/**
 * Destructor for HMM. Sets previously allocated memory of all class attributes free.
 *
 */
        virtual ~InitialProbability();

        double* getInitialProb();
        void updateSample(double** gamma, int i);
        void updateSampleCoupled(double** gamma, int i, int* couples, SEXP bidirOptimParams, int* T, int currN);
        void update(int nsample, SEXP bidirOptimParams, int* T);
        int getK();
        void finalize();
};
#endif
