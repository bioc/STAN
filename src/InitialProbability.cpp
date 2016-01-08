
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


void InitialProbability::updateSampleCoupled(double** gamma, int i, int* couples, SEXP bidirOptimParams, int* T, int currN)
{
    if(LENGTH(bidirOptimParams) == 0)
    {
//printf("here->initial\n");
        int twin_for_i = couples[i];
        int t;
        for(t=1; t<T[currN]; t++)
        {
            this->updateNumeratorPI[i] += (gamma[t-1][i]+gamma[t][twin_for_i]);
        }

//	printf("%d: %f\n", i, this->updateNumeratorPI[i]);
    }
    else if(LENGTH(bidirOptimParams) != 0)        // numerical optimization of initial state probabilities
    {
        this->updateNumeratorPI[i] +=  gamma[0][i];
        REAL(getListElement(bidirOptimParams, "initGamma"))[i] += gamma[0][i];
//Rprintf("gamma[%d]=%f\n",i,REAL(getListElement(bidirOptimParams, "initGamma"))[i]);
    }
/*else {
    this->updateNumeratorPI[i] += gamma[0][i];
}*/
//this->finalize();
}


void InitialProbability::update(int nsample, SEXP bidirOptimParams, int* T)
{
    int i;
    for(i=0; i<this->K; i++)
    {
        if(LENGTH(bidirOptimParams) != 0)
        {
            if(INTEGER(getListElement(bidirOptimParams, "update"))[0] == 1)
            {
//Rprintf("%d ", INTEGER(getListElement(bidirOptimParams, "update"))[0]);
                this->pi[i] = REAL(getListElement(bidirOptimParams, "statD"))[i];
            }
            REAL(getListElement(bidirOptimParams, "initGamma"))[i] = 0;
            this->updateNumeratorPI[i] = 0;
//	Rprintf("%f ", this->pi[i]);
        }
        else if(T != NULL)
        {
            int allT = 0;
            int n;
            for(n=0; n<nsample; n++)
            {
                allT += T[n];
            }
            this->pi[i] = this->updateNumeratorPI[i]/(2*allT-2);
            this->updateNumeratorPI[i] = 0;
        }
        else
        {
//printf("here-init->update\n");
            this->pi[i] = this->updateNumeratorPI[i]/nsample;
            this->updateNumeratorPI[i] = 0;
        }
//Rprintf("\n");
    }
//this->finalize();

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
