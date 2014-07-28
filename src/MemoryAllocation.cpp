
#include "MemoryAllocation.h"
#include "ParamContainerEmissions.h"

double** allocateNumericMatrix(int d1, int d2)
{
    double **matrix = (double**)malloc(sizeof(double*)*d1);
    if(matrix == NULL)
    {
        error("Not enough memory!\n");
    }

    int i;
    for(i=0; i<d1; i++)
    {
        matrix[i] = allocateNumericVector(d2);
    }
    return matrix;
}


double* allocateNumericVector(int d)
{
    double *vector = (double*)malloc(sizeof(double)*d);
    if(vector == NULL)
    {
        error("Not enough memory!\n");
    }
    return vector;
}


int allocateMemAlpha(double*** alpha, int maxLen, int K)
{
    int t,i;
    int memory_used = 0;

    *alpha = (double**)malloc(sizeof(double*)*maxLen);
    memory_used += sizeof(double*)*maxLen;

    for(t=0; t<maxLen; t++)
    {
        (*alpha)[t] = (double*)malloc(sizeof(double)*K);
        memory_used += sizeof(double)*K;
        for(i=0; i<K; i++)
        {
            (*alpha)[t][i] = 0;
        }
    }
    double megabytes_used = ((double)memory_used)/1000000;

    if(DEBUG_MEMORY)
    {
        printf("Alpha needs %lf MB of memory.\n", megabytes_used);
    }
    return memory_used;
}


int allocateMemBeta(double*** beta, int maxLen, int K)
{
    int t,i;
    int memory_used = 0;

    *beta  = (double**)malloc(sizeof(double*)*maxLen);
    memory_used += sizeof(double*)*maxLen;

    for(t=0; t<maxLen; t++)
    {
        (*beta)[t] = (double*)malloc(sizeof(double)*K);
        for(i=0; i<K; i++)
        {
            (*beta)[t][i] = 0;
        }
    }

    double megabytes_used = ((double)memory_used)/1000000;

    if(DEBUG_MEMORY)
    {
        printf("Beta needs %lf MB of memory.\n", megabytes_used);
    }
    return memory_used;

}


int allocateMemRescFac(double** c, int maxLen, int K)
{
    int i,t;

    int memory_used = 0;

    *c = (double*)malloc(sizeof(double)*maxLen);
    memory_used += sizeof(double)*maxLen;
    for(t=0; t<maxLen; t++)
    {
        (*c)[t] = 0;
    }
    double megabytes_used = ((double)memory_used)/1000000;

    if(DEBUG_MEMORY)
    {
        printf("Rescaling factor needs %lf MB of memory.\n", megabytes_used);
    }
    return memory_used;

}


int allocateMemGamma(double*** gamma, int maxLen, int K)
{
    int i,t;

    int memory_used = 0;

    *gamma  = (double**)malloc(sizeof(double*)*maxLen);
    memory_used += sizeof(double*)*maxLen;

    for(t=0; t<maxLen; t++)
    {
        (*gamma)[t] = (double*)malloc(sizeof(double)*K);
        memory_used += sizeof(double)*K;
        for(i=0; i<K; i++)
        {
            (*gamma)[t][i] = 0;
        }
    }

    double megabytes_used = ((double)memory_used)/1000000;

    if(DEBUG_MEMORY)
    {
        printf("Gamma needs %lf MB of memory.\n", megabytes_used);
    }
    return memory_used;

}


int allocateMemXsi(double**** xsi, int maxLen, int K)
{
    int i,j,t;

    int memory_used = 0;

    *xsi  = (double***)malloc(sizeof(double**)*maxLen);
    memory_used += sizeof(double**)*maxLen;

    for(t=0; t<maxLen; t++)
    {
        (*xsi)[t] = (double**)malloc(sizeof(double*)*K);
        memory_used += sizeof(double*)*K;

        for(i=0; i<K; i++)
        {
            (*xsi)[t][i] = (double*)malloc(sizeof(double)*K);
            memory_used += sizeof(double)*K;
            for(j=0; j<K; j++)
            {
                (*xsi)[t][i][j] = 0;
            }
        }
    }

    double megabytes_used = ((double)memory_used)/1000000;

    if(DEBUG_MEMORY)
    {
        printf("Xsi needs %lf MB of memory.\n", megabytes_used);
    }
    return memory_used;

}


int allocateMemEmissionProb(double*** emissionProb, int maxLen, int K)
{
    int i,j,t;

    int memory_used = 0;

    *emissionProb = (double**)malloc(sizeof(double*)*K);
    memory_used += sizeof(double*)*K;

    for(i=0; i<K; i++)
    {
        (*emissionProb)[i] = (double*)malloc(sizeof(double)*maxLen);
        memory_used += sizeof(double)*maxLen;
        for(t=0; t<maxLen; t++)
        {
            (*emissionProb)[i][t] = 0;
        }
    }

    double megabytes_used = ((double)memory_used)/1000000;

    if(DEBUG_MEMORY)
    {
        printf("Emission probabilities need %lf MB of memory.\n", megabytes_used);
    }
    return memory_used;

}


int deallocateMemAlpha(double** alpha, int maxLen, int K)
{
    int i,t;
    int memory_free = 0;

    for(t=0; t<maxLen; t++)
    {
        free(alpha[t]);
        memory_free += sizeof(double)*K;
    }

    free(alpha);
    memory_free += sizeof(double*)*maxLen;

    return memory_free;
}


int deallocateMemBeta(double** beta, int maxLen, int K)
{
    int i,j,t;
    int memory_free = 0;

    for(t=0; t<maxLen; t++)
    {
        free(beta[t]);
        memory_free += sizeof(double)*K;
    }

    free(beta);
    memory_free += sizeof(double*)*maxLen;

    return memory_free;

}


int deallocateMemRescFac(double* c, int maxLen, int K)
{
    int i,j,t;
    int memory_free = 0;

    free(c);
    memory_free += sizeof(double)*maxLen;

    return memory_free;

}


int deallocateMemGamma(double** gamma, int maxLen, int K)
{
    int i,j,t;
    int memory_free = 0;

    for(t=0; t<maxLen; t++)
    {
        free(gamma[t]);
        memory_free += sizeof(double)*K;

    }

    free(gamma);
    memory_free += sizeof(double*)*maxLen;

    return memory_free;

}


int deallocateMemXsi(double*** xsi, int maxLen, int K)
{
    int i,j,t;
    int memory_free = 0;

    for(t=0; t<maxLen; t++)
    {
        for(i=0; i<K; i++)
        {
            free(xsi[t][i]);
            memory_free += sizeof(double)*K;
        }
        free(xsi[t]);
        memory_free += sizeof(double)*K;
    }

    free(xsi);
    memory_free += sizeof(double**)*maxLen;

    return memory_free;

}


int deallocateMemEmissionProb(double** emissionProb, int maxLen, int K)
{
    int i,j,t;
    int memory_free = 0;

    for(i=0; i<K; i++)
    {
        free(emissionProb[i]);
        memory_free += sizeof(double)*maxLen;
    }
    free(emissionProb);

    return memory_free;

}
