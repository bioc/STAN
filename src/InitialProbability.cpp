
#include "InitialProbability.h"

InitialProbability::InitialProbability(double *pi, int K)
{
    this->pi = pi;
    this->K = K;
    this->updateNumeratorPI = (double*)malloc(sizeof(double)*K);

    int i;
    for(i=0; i<K; i++)
    {
        this->updateNumeratorPI[i] = 0;
    }
    if(DEBUG_MEMORY)
    {
        printf("new->InitialProbability; (%d bytes)\n", sizeof(double)*K*2);
    }

}


InitialProbability::~InitialProbability()
{
    free(this->pi);
    free(this->updateNumeratorPI);

    if(DEBUG_MEMORY)
    {
        printf("delete->InitialProbability; (%d bytes)\n", sizeof(double)*K*2);
    }
}


double* InitialProbability::getInitialProb()
{
    return this->pi;
}


void InitialProbability::updateSample(double** gamma, int i)
{
    this->updateNumeratorPI[i] += gamma[0][i];
}


void InitialProbability::updateSampleCoupled(double** gamma, int i, int* couples, SEXP bidirOptimParams)
{
    if(LENGTH(bidirOptimParams) == 0)
    {
        int twin_for_i = couples[i];

        this->updateNumeratorPI[i] += (gamma[0][i]+gamma[0][twin_for_i])/2;
    }
    // numerical optimization of initial state probabilities
    else if(LENGTH(bidirOptimParams) != 0)       
    {
        this->updateNumeratorPI[i] +=  gamma[0][i];
        REAL(getListElement(bidirOptimParams, "initGamma"))[i] += gamma[0][i];
    }
}


void InitialProbability::update(int nsample, SEXP bidirOptimParams)
{
    int i;
    for(i=0; i<this->K; i++)
    {
        if(LENGTH(bidirOptimParams) != 0)
        {
            if(INTEGER(getListElement(bidirOptimParams, "update"))[0] == 1)
            {
                this->pi[i] = REAL(getListElement(bidirOptimParams, "statD"))[i];
            }
            REAL(getListElement(bidirOptimParams, "initGamma"))[i] = 0;
            this->updateNumeratorPI[i] = 0;
        }
        else
        {
            this->pi[i] = this->updateNumeratorPI[i]/nsample;
            this->updateNumeratorPI[i] = 0;
        }
    }
}


int InitialProbability::getK()
{
    return this->K;
}


void InitialProbability::finalize()
{
    double sum=0;
    int i;
    for(i=0; i<this->K; i++)
    {
        sum += this->pi[i];
    }

    for(i=0; i<this->K; i++)
    {
        this->pi[i] /= sum;
    }
}
