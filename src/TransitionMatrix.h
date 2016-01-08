
#ifndef TRANSMAT_HEADER
#define TRANSMAT_HEADER

#include <stdio.h>
#include <stdlib.h>
#include "DebugConstants.h"
#include <R.h>
#include <Rdefines.h>
#include "RAccessUtils.h"

using namespace std;

class TransitionMatrix
{
    protected:

        int K;
        double** A;                               //!<@brief KxK matrix of transition probabilities between states.
        double** updateNumerator;
        double** updateDenominator;
// TODO: also constraints here!

    public:

        TransitionMatrix(double **A, int K);

/**
 * Destructor for HMM. Sets previously allocated memory of all class attributes free.
 *
 */
        virtual ~TransitionMatrix();

        double** getTransMat();
        int getK();
        void update(double effective_zero);
        void update(int* couples, double effective_zero);
        void update(SEXP bidirOptimParams);
        void updateAuxiliaries(double** gamma, double*** xsi,  double* Pk, int* T, int n, int** isNaN, int ncores, double effective_zero, int verbose);
        void updateAuxiliaries(double** gamma, double*** xsi,  double* Pk, int* T, int n, int* couples, SEXP bidirOptimParams, int** isNaN, int ncores, double effective_zero, int verbose);
        void finalize();
        SEXP callRsolnp(SEXP pars);
};
#endif
