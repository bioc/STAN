
/**
 * @file    HMM.h
 * @author  Benedikt Zacher, AG Tresch, Gene Center Munich (zacher@lmb.uni-muenchen.de)
 * @version 0.99.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 *
 * TODO
 */

#ifndef HMM_HEADER
#define HMM_HEADER

#define CFree(x) if(x != NULL) free(x);

#include <new>
#include <list>
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <R_ext/Utils.h>
#include "EmissionFunction.h"
#include "TransitionMatrix.h"
#include "InitialProbability.h"
#include "DebugConstants.h"

using namespace std;

class HMM
{
protected:

    int K; 
    InitialProbability* pi;
    TransitionMatrix* A;
    EmissionFunction** emissions;

public:

    HMM() {}
    HMM(int K, InitialProbability* pi, TransitionMatrix* A, EmissionFunction** emissions);
    ~HMM();  
    
    virtual void getAlphaBeta(double*** obs, double **alpha, double **beta, double *c, double **emissionProb, int* T, int n, int ncores, double effective_zero, int verbose);
    list<double> BaumWelch(double*** observations, int* T, int nsample, int maxIters, int** flags, int* state2flag, int* couples, int* revop, int verbose, int updateTransMat, int** isNaN, double*** fixedEmission, SEXP bidirOptimParams, SEXP emissionPrior, int ncores, double effective_zero, double convergence, int incrementalEM);
    virtual int allocateMemEM(double*** emissionProb, double*** alpha, double*** beta, double*** gamma, double**** xsi, double** c, double** Pk, int maxLen, int nsample);
    virtual int deallocateMemEM(double** emissionProb, double** alpha, double** beta, double** gamma, double*** xsi, double* c, double* Pk, int maxLen, int nsample);
    virtual void getGammaXsi(double** alpha, double** beta, double* c, double** emissionProb, double** gamma, double*** xsi, int* T, int n, int ncores, double effective_zero, int verbose);
    virtual void calcEmissionProbs(double*** obs, double** emissionProb, int* T, int n, int** flags, int* state2flag, int* revop, int** isNaN, int ncores, int verbose);
    void reverseObs(double *orig, double** rev, int* revop, int D);
    void Viterbi(int **S, double*** obs, int nsample, int* T, int verbose, int** isNaN, double*** fixedEmission);
    void updateSampleAux(double*** observations, int* T, int n, double** alpha, double** beta, double** gamma, double*** xsi, double* Pk, int* state2flag, int* couples, int* revop, int** isNaN, double*** fixedEmission, SEXP bidirOptimParams, SEXP emissionPrior, int ncores, double effective_zero, int verbose);
    void getGamma(double** alpha, double** beta, double* c,  double** emissionProb, double** gamma, int* T, int n, int ncores, double effective_zero, int verbose);
    void updateModelParams(double*** observations, int nsample, int* state2flag, int* couples, int* revop, int verbose, int updateTransMat, int** isNaN, double*** fixedEmission, SEXP bidirOptimParams, SEXP emissionPrior, int ncores, double effective_zero, int* myStateBuckets, double* Pk, int curriter, int currN);

};
#endif
