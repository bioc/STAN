
#include "ParamContainerEmissions.h"
#define ARRAYSIZE(a) (sizeof(a) / sizeof(a[0]))

int ParamContainerEmissions::getWhichOne()
{
    return this->whichone;
}


ParamContainerEmissions** allocateParamContainerVector(int d)
{
    ParamContainerEmissions **vector = (ParamContainerEmissions**)malloc(sizeof(ParamContainerEmissions*)*d);
    if(vector == NULL)
    {
        error("Not enough memory!\n");
    }
    return vector;
}


ParamContainerEmissions::ParamContainerEmissions(int d)
{
    this->D=d;
    this->whichone = JOINTLYINDEPENDENT;

}


//ParamContainerEmissions::~ParamContainerEmissions() {}

/**
 * 01) ParamContainerEmissions Constructor for MULTIVARIATEGAUSSIAN
 *
 */
ParamContainerEmissions::ParamContainerEmissions(double **mu, double **sigma, double regularize, int D, int* start, int updateCov, int sharedCov)
{
    this->parallelUpdate = 1;
    this->logCovPrior = 0;
    this->updateCov = updateCov;
    this->sharedCov = sharedCov;
    this->whichone = MULTIVARIATEGAUSSIAN;
    this->mu = mu;
    this->sigma = sigma;
    this->regularize = regularize;
    this->D = D;
    this->start=start;
    this->inverseSigma = allocateNumericMatrix(D, D);
    int i,j;
    for(i=0; i<D; i++)
    {
        for(j=0; j<D; j++)
        {
            this->inverseSigma[i][j] = this->sigma[i][j];
//printf("%f ", this->inverseSigma[i][j]);
        }
//printf("\n");
    }
//printf("\n");
    inverse(this->inverseSigma, D);
//printf("whichone: %d dimension %d , arraystartsize %d \n", this->whichone, D, ARRAYSIZE(start));
//LapackInvAndDet(this->inverseSigma, D);

    this->determinant = matrixDet(sigma, D);
//printf("determinant %d ", this->determinant);
    int mem = sizeof(double*)*D + sizeof(double)*1 + 2*sizeof(double*)*D + 2*sizeof(double)*D;
    if(DEBUG_MEMORY)
    {
        printf("new->ParamContainerEmissions:MULTIVARIATEGAUSSIAN ; (%d bytes) ", mem);
    }

}


/**
 * 02) ParamContainerEmissions Constructor for BERNOULLI
 *
 */

ParamContainerEmissions::ParamContainerEmissions(double p, int D, int* start)
{
    this->parallelUpdate = 1;

    this->p = p;
    this->whichone = BERNOULLI;
    this->D = D;
    this->start=start;
    int mem = 2*sizeof(double)*D;

    if(DEBUG_MEMORY)
    {
        printf("new->ParamContainerEmissions:BERNOULLI; (%d bytes) ", mem);
    }
}


ParamContainerEmissions::ParamContainerEmissions(double lambda, int D, int* start, int whichone)
{
    this->parallelUpdate = 1;

    if(whichone != POISSON)
    {
        error("Must be Poisson emission here!\n");
    }
    this->lambda = lambda;
    this->whichone = POISSON;
    this->D = D;
    this->start=start;
    int mem = 2*sizeof(double)*D;

    if(DEBUG_MEMORY)
    {
        printf("new->ParamContainerEmissions:POISSON; (%d bytes) ", mem);
    }
}


ParamContainerEmissions::ParamContainerEmissions(double* p, int* reverseComplementary, int D, int* start, int stateDir)
{
    this->parallelUpdate = 1;
    this->mp = p;
    this->reverseComplementary = reverseComplementary;
    this->whichone = MULTINOMIAL;
    this->D = D;
    this->stateDir = -stateDir;
//Rprintf("stateDir=%d\n", this->stateDir);
    this->start=start;
    int mem = 2*sizeof(double)*D;

    if(DEBUG_MEMORY)
    {
        printf("new->ParamContainerEmissions:MULTINOMIAL; (%d bytes) ", mem);
    }
}


ParamContainerEmissions::ParamContainerEmissions(double mu_nb, double size_nb, double* sizeFactor_nb, double pi_nb, int D, int* start, SEXP uniqueCountSplit)
{

    this->parallelUpdate = 0;
    this->mu_nb = mu_nb;
    this->size_nb = size_nb;
    this->sizeFactor_nb = sizeFactor_nb;
    this->pi_nb = pi_nb;
    this->whichone = NEGATIVEBINOMIAL;
    this->D = D;
    this->start=start;
    this->uniqueCountSplit = uniqueCountSplit;
//Rprintf("mu=%f, size=%f\n", this->mu_nb, this->size_nb);
    int mem = 2*sizeof(double)*D;

    if(DEBUG_MEMORY)
    {
        printf("new->ParamContainerEmissions:NEGATIVEBINOMIAL; (%d bytes) ", mem);
    }
}


ParamContainerEmissions::ParamContainerEmissions(double mu_poilog, double sigma_poilog, double* sigma_poilog_sf, double* sizeFactor_poilog, int D, int* start, SEXP uniqueCountSplit)
{

    this->parallelUpdate = 0;
    this->mu_poilog = mu_poilog;
    this->sigma_poilog = sigma_poilog;
    this->sigma_poilog_sf = sigma_poilog_sf;
    this->sizeFactor_poilog = sizeFactor_poilog;
    this->whichone = POISSONLOGNORMAL;
    this->D = D;
    this->start=start;
    this->uniqueCountSplit = uniqueCountSplit;
//Rprintf("mu=%f, size=%f\n", this->mu_nb, this->size_nb);
    int mem = 2*sizeof(double)*D;

    if(DEBUG_MEMORY)
    {
        printf("new->ParamContainerEmissions:POISSONLOGNORMAL; (%d bytes) ", mem);
    }
}


/**
 * ParamContainerEmissions DESTRUCTOR
 */
ParamContainerEmissions::~ParamContainerEmissions()
{
//printf("I AM HERE!!!\n");
    int mem = sizeof(double*)*this->D;

    if(this->whichone == MULTIVARIATEGAUSSIAN)
    {
        int i;
        for(i=0; i<this->D; i++)
        {
            free(this->mu[i]);
            free(this->sigma[i]);
            free(this->inverseSigma[i]);
        }
        free(this->mu);
        free(this->sigma);
        free(this->inverseSigma);

        mem = sizeof(double*)*this->D + sizeof(double)*1 + 2*sizeof(double*)*this->D + 2*sizeof(double)*this->D;

        if(DEBUG_MEMORY)
        {
            printf("delete->ParamContainerEmissions:MULTIVARIATEGAUSSIAN; (%d bytes) \n", mem);
        }

    }
    if(this->whichone == BERNOULLI)
    {
        int i;

        if(DEBUG_MEMORY)
        {
            printf("delete->ParamContainerEmissions:BERNOULLI; (%d bytes) \n", mem);
        }
    }
    if(this->whichone == POISSON)
    {
        int i;

        if(DEBUG_MEMORY)
        {
            printf("delete->ParamContainerEmissions:BERNOULLI; (%d bytes) \n", mem);
        }
    }
    if(this->whichone == MULTINOMIAL)
    {
        if(DEBUG_MEMORY)
        {
            printf("delete->ParamContainerEmissions:MULTINOMIAL; (%d bytes) \n", mem);
        }
        free(mp);
        free(reverseComplementary);
    }
    if(this->whichone == NEGATIVEBINOMIAL)
    {
        if(DEBUG_MEMORY)
        {
            printf("delete->ParamContainerEmissions:NEGATIVEBINOMIAL; (%d bytes) \n", mem);
        }

        free(this->sizeFactor_nb);

    }
    if(this->whichone == POISSONLOGNORMAL)
    {
        if(DEBUG_MEMORY)
        {
            printf("delete->ParamContainerEmissions:POISSONLOGNORMAL; (%d bytes) \n", mem);
        }

        free(this->sizeFactor_poilog);

    }

    if(this->whichone == POISSONLOGNORMAL || this->whichone == NEGATIVEBINOMIAL)
    {
        int n;
        for(n=0; n<this->nsample; n++)
        {
            free(this->uniqueObsProb[n]);
            free(this->myUniqueLens[n]);
        }
        free(this->uniqueObsProb);
        free(this->myUniqueLens);
//		free(start);

    }

    if(this->whichone == JOINTLYINDEPENDENT)
    {

        int n;
        for(n=0; n<this->nsample; n++)
        {
            free(this->gammaAux[n]);
        }
        free(this->gammaAux);
//free(T);

    }

}


/*
 * GENERAL GETTER and SETTER
 */
void ParamContainerEmissions::setDataVars(int nsample, int* T)
{
    int n,t;

    this->nsample = nsample;
    this->T = T;
    this->gammaAux = (double**)malloc(sizeof(double*)*nsample);
    for(n=0; n<nsample; n++)
    {
        this->gammaAux[n] = (double*)malloc(sizeof(double)*T[n]);
        for(t=0; t<T[n]; t++)
        {
            this->gammaAux[n][t] = 0;
        }
    }
}


void ParamContainerEmissions::initUniqueObsProb(double*** observations, int nsample, int* T, int* revop)
{
    int t, d, n, i;
    double myMax;

    this->myUniqueLens = (int**)malloc(sizeof(int*)*nsample);
    this->uniqueObsProb = (double**)malloc(sizeof(double*)*nsample);

    for(n=0; n<nsample; n++)
    {
        this->myUniqueLens[n] = (int*)malloc(sizeof(int)*this->D);
        for(d=0; d<this->D; d++)
        {
            int otherDim = this->start[d];
            if(revop != NULL)
            {
                otherDim = revop[this->start[d]];
            }
            myMax=0;
            for(t=0; t<T[n]; t++)
            {
                if(observations[n][t][this->start[d]] == observations[n][t][this->start[d]])
                {
                    if(myMax < observations[n][t][this->start[d]])
                    {
                        myMax = observations[n][t][this->start[d]];
                    }
                    if(myMax < observations[n][t][otherDim])
                    {
                        myMax = observations[n][t][otherDim];
                    }
                }
            }
            myMax = myMax+1;
            this->myUniqueLens[n][d] = (int)myMax;
//	Rprintf("n=%d, d=%d, #unique=%d\n", n, this->start[d], this->myUniqueLens[n][d]);
            this->uniqueObsProb[n] = (double*)malloc(sizeof(double)*myMax);
            for(i=0; i<this->myUniqueLens[n][d]; i++)
            {
                this->uniqueObsProb[n][i] = -1;
            }
            for(t=0; t<T[n]; t++)
            {
                if(observations[n][t][this->start[d]] == observations[n][t][this->start[d]])
                {
                    int currIndex = (int)(observations[n][t][this->start[d]]);
                    this->uniqueObsProb[n][currIndex] = 1;
                    int otherDim = this->start[d];
                    if(revop != NULL)
                    {
                        otherDim = revop[this->start[d]];
                        currIndex = (int)(observations[n][t][otherDim]);
                        this->uniqueObsProb[n][currIndex] = 1;
                    }
                }
            }
        }
    }

}


double** ParamContainerEmissions::getUniqueObsProb()
{
    return this->uniqueObsProb;
}


int** ParamContainerEmissions::getUniqueLens()
{
    return this->myUniqueLens;
}


/*
 * GENERAL GETTER and SETTER
 */
void ParamContainerEmissions::setDataVars(double** wrapper_gamma, int nsample, int* T)
{
//int n,t;

    this->nsample = nsample;
    this->T = T;
    this->gammaAux = wrapper_gamma;
//int n;
/*for(n=0; n<this->nsample; n++) {
    this->gammaAux[n] = wrapper_gamma[n];
}*/
}


int* ParamContainerEmissions::getT()
{
    return this->T;
}


int ParamContainerEmissions::getNsample()
{
    return this->nsample;
}


int ParamContainerEmissions::getCurrState()
{
    return this->currState;
}


void ParamContainerEmissions::setCurrState(int state)
{
    this->currState = state;
}


int ParamContainerEmissions::getD()
{
    return this->D;
}


int* ParamContainerEmissions::getStart()
{
    return this->start;
}


/*
 * MULTIVARIATEGAUSSIAN GETTER and SETTER
 */

void ParamContainerEmissions::setGammaAux(double val, int n, int t)
{
    this->gammaAux[n][t] = val;
}


double** ParamContainerEmissions::getGammaAux()
{
    return this->gammaAux;
}


double** ParamContainerEmissions::getGaussianMU()
{
    return this->mu;
}


double** ParamContainerEmissions::getGaussianSIGMA()
{
    return this->sigma;
}


double** ParamContainerEmissions::getGaussianINVSIGMA()
{
    return this->inverseSigma;
}


double ParamContainerEmissions::getGaussianDET()
{
    return this->determinant;
}


double ParamContainerEmissions::getGaussianREG()
{
    return this->regularize;
}


int ParamContainerEmissions::getUpdateCov()
{
    return this->updateCov;
}


void ParamContainerEmissions::setLogCovPrior(double prior)
{
    this->logCovPrior = prior;
}


double ParamContainerEmissions::getLogCovPrior()
{
    return this->logCovPrior;
}


void ParamContainerEmissions::setGaussianMU(double **mu)
{
    int i;
    for(i=0; i<D; i++)
    {
        this->mu[i][0] = mu[i][0];
    }
}


void ParamContainerEmissions::setGaussianMUelement(double val, int d)
{
    this->mu[d][0] = val;
}


void ParamContainerEmissions::setGaussianSIGMAelement(double val, int d1, int d2)
{
    this->sigma[d1][d2] = val;
}


void ParamContainerEmissions::setGaussianINVSIGMAelement(double val, int d1, int d2)
{
    this->inverseSigma[d1][d2] = val;
}


void ParamContainerEmissions::setGaussianDET(double val)
{
    this->determinant = val;
}


void ParamContainerEmissions::setGaussianSIGMA(double **sigma)
{
    int i,j;
    for(i=0; i<D; i++)
    {
        for(j=0; j<D; j++)
        {
            this->sigma[i][j] = sigma[i][j];
            this->inverseSigma[i][j] = sigma[i][j];
        }
    }
    inverse(inverseSigma, this->D);
    this->determinant = matrixDet(sigma, this->D);
}


/*
 * BERNOULLI GETTER and SETTER
 *
 */

double ParamContainerEmissions::getBernoulliP()
{
    return this->p;
}


void ParamContainerEmissions::setBernoulliPelement(double val, int d)
{
    this->p = val;

}


void ParamContainerEmissions::setBernoulliP(double p)
{
    this->p = p;

}


double  ParamContainerEmissions::getPoissonLambda()
{
    return this->lambda;
}


void  ParamContainerEmissions::setPoissonLambdaelement(double val, int d)
{
    this->lambda = val;
}


void  ParamContainerEmissions::setPoissonLambda(double lambda)
{
    this->lambda = lambda;
}


double* ParamContainerEmissions::getMultinomialP()
{
    return this->mp;
}


int ParamContainerEmissions::getN()
{
    return this->n;
}


int* ParamContainerEmissions::getReverseComplementary()
{
    return this->reverseComplementary;
}


int ParamContainerEmissions::getStateDir()
{
    return this->stateDir;
}


void ParamContainerEmissions::setMultinomialPelement(double val, int d)
{
    this->mp[d] = val;
}


void ParamContainerEmissions::setMultinomialP(double* p)
{
    this->mp = p;
}


double ParamContainerEmissions::getMuNB()
{
    return this->mu_nb;
}


double ParamContainerEmissions::getSizeNB()
{
    return this->size_nb;
}


double ParamContainerEmissions::getPiNB()
{
    return this->pi_nb;
}


double ParamContainerEmissions::setMuNB(double mu_nb)
{
    this->mu_nb = mu_nb;
}


double ParamContainerEmissions::setSizeNB(double size_nb)
{
    this->size_nb = size_nb;
}


double ParamContainerEmissions::setPiNB(double pi_nb)
{
    this->pi_nb = pi_nb;
}


double* ParamContainerEmissions::getSizeFactorNB()
{
    return this->sizeFactor_nb;
}


int ParamContainerEmissions::getSharedCov()
{
    return this->sharedCov;
}


double ParamContainerEmissions::getMuPoiLog()
{
    return this->mu_poilog;
}


double ParamContainerEmissions::getSigmaPoiLog()
{
    return this->sigma_poilog;
}


double ParamContainerEmissions::setMuPoiLog(double mu_poilog)
{
    this->mu_poilog = mu_poilog;
}


double ParamContainerEmissions::setSigmaPoiLog(double sigma_poilog)
{
    this->sigma_poilog = sigma_poilog;
}


double ParamContainerEmissions::getSigmaPoiLogVec(int n)
{
    return this->sigma_poilog_sf[n];
}


double* ParamContainerEmissions::getSigmaPoiLogVec()
{
    return this->sigma_poilog_sf;
}


void ParamContainerEmissions::setSigmaPoiLogVec(double val, int n)
{
    this->sigma_poilog_sf[n] = val;
}


double* ParamContainerEmissions::getSizeFactorPoiLog()
{
    return this->sizeFactor_poilog;
}


SEXP ParamContainerEmissions::getUniqueCountSplit()
{
    return this->uniqueCountSplit;
}


int ParamContainerEmissions::doParallelUpdate()
{
    return this->parallelUpdate;
}
