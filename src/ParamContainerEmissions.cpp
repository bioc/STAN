
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


ParamContainerEmissions::ParamContainerEmissions(double **mu, double **sigma, double regularize, int D, int* start)
{
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
        }
    }
    inverse(this->inverseSigma, D);
    this->determinant = matrixDet(sigma, D);
    int mem = sizeof(double*)*D + sizeof(double)*1 + 2*sizeof(double*)*D + 2*sizeof(double)*D;
    if(DEBUG_MEMORY)
    {
        printf("new->ParamContainerEmissions:MULTIVARIATEGAUSSIAN ; (%d bytes) ", mem);
    }

}


ParamContainerEmissions::ParamContainerEmissions(double p, int D, int* start)
{

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


ParamContainerEmissions::~ParamContainerEmissions()
{
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

        int mem = sizeof(double*)*this->D + sizeof(double)*1 + 2*sizeof(double*)*this->D + 2*sizeof(double)*this->D;

        if(DEBUG_MEMORY)
        {
            printf("delete->ParamContainerEmissions:MULTIVARIATEGAUSSIAN; (%d bytes) \n", mem);
        }

    }
    if(this->whichone == BERNOULLI)
    {
        int i;

        int mem = sizeof(double*)*this->D;
        if(DEBUG_MEMORY)
        {
            printf("delete->ParamContainerEmissions:BERNOULLI; (%d bytes) \n", mem);
        }
    }

}


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


void ParamContainerEmissions::setDataVars(double** wrapper_gamma)
{
    this->gammaAux = wrapper_gamma;
}


int* ParamContainerEmissions::getT()
{
    return this->T;
}


int ParamContainerEmissions::getNsample()
{
    return this->nsample;
}


int ParamContainerEmissions::getD()
{
    return this->D;
}


int* ParamContainerEmissions::getStart()
{
    return this->start;
}


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
