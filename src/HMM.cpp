
#include "HMM.h"

HMM::HMM(int K, InitialProbability* pi, TransitionMatrix* A, EmissionFunction** emissions)
{
    if(DEBUG_MEMORY)
    {
        Rprintf("new->HMM;\n");
    }
    this->K = K;
    this->pi = pi;
    this->A = A;
    this->emissions = emissions;

}


HMM::~HMM()
{
    if(DEBUG_MEMORY)
    {
        Rprintf("delete->HMM;\n");
    }
    delete this->pi;
    delete this->A;
    int i;

    if(this->emissions != NULL)
    {
        for(i=0; i<K; i++)
        {
            if(DEBUG_MEMORY)
            {
                Rprintf("\tstate=%d: ", i);
            }
            delete this->emissions[i];
        }

        free(this->emissions);
    }

}


void HMM::reverseObs(double *orig, double** rev, int* revop, int D)
{
    int i;
    for(i=0; i<D; i++)
    {
        (*rev)[i] = orig[revop[i]];
    }
}


void HMM::calcEmissionProbs(double*** obs, double** emissionProb, int* T, int n, int** flags, int* state2flag, int* revop, int** isNaN, int ncores, int verbose)
{
    double proba;
    int D = this->emissions[0]->getParameter()->getD();
    double* rev_obs = (double*)malloc(sizeof(double)*D);

    if(verbose)
    {
        Rprintf("Sequence %d => Calculating emission probabilities.                                                                           \r", n+1);
    }

    if(DEBUG)
    {
        // flags, strand-specific observations
        if(state2flag != NULL && flags != NULL && revop != NULL)
        {
            Rprintf("density calculation: +flags +reverse observations\n");
        }
        // strand-specific observation, no flags
        else if(state2flag != NULL && flags == NULL && revop != NULL)
        {
            Rprintf("density calculation: -flags +reverse observations\n");
        }
        else if(state2flag != NULL && flags != NULL && revop == NULL)
        {
            Rprintf("density calculation: +flags -reverse observations\n");
        }
        // no flags, no strand-specific observations
        else                                     
        {
            Rprintf("density calculation: -flags -reverse observations\n");
        }
    }

    int seglen = (int) T[n]/ncores;
    int* myBounds = (int*)malloc(sizeof(int)*(ncores+1));

    myBounds[0] = 0;
    int iter;
    for(iter=1; iter<ncores; iter++)
    {
        myBounds[iter] = myBounds[iter-1]+seglen;
    }
    myBounds[ncores] = T[n];

    int curr_core;
    double **myTransMat = this->A->getTransMat();

    #pragma omp parallel for num_threads(ncores)
    for(curr_core=0; curr_core<ncores; curr_core++)
    {
        int i,j,t,k;
        double denom;
        for(t=myBounds[curr_core]; t<myBounds[curr_core+1]; t++)
        {
            for(i=0; i<this->K; i++)
            {
                // flags, strand-specific observations
                if(state2flag != NULL && flags != NULL && revop != NULL)
                {
                    if(state2flag[i] == 1 && (flags[n][t] != 1))
                    {
                        reverseObs(obs[n][t], &rev_obs, revop, D);
                        emissionProb[i][t] = this->emissions[i]->calcEmissionProbability(rev_obs, isNaN[n][t]);
                    }
                    else if(state2flag[i] == -1 && (flags[n][t] != -1))
                    {
                        emissionProb[i][t] = this->emissions[i]->calcEmissionProbability(obs[n][t], isNaN[n][t]);
                    }
                    else if(state2flag[i] == -100)
                    {
                        emissionProb[i][t] = this->emissions[i]->calcEmissionProbability(obs[n][t], isNaN[n][t]);
                    }
                    else
                    {
                        emissionProb[i][t] = 0;
                    }
                }
                // strand-specific observation, no flags
                else if(state2flag != NULL && flags == NULL && revop != NULL)
                {
                    if(state2flag[i] == 1)
                    {
                        reverseObs(obs[n][t], &rev_obs, revop, D);
                        emissionProb[i][t] = this->emissions[i]->calcEmissionProbability(rev_obs, isNaN[n][t]);
                    }
                    else
                    {
                        emissionProb[i][t] = this->emissions[i]->calcEmissionProbability(obs[n][t], isNaN[n][t]);
                    }

                }
                // flags, no strand-specific observations
                else if(state2flag != NULL && flags != NULL && revop == NULL)
                {
                    if(state2flag[i] != flags[n][t])
                    {
                        emissionProb[i][t] = this->emissions[i]->calcEmissionProbability(obs[n][t], isNaN[n][t]);
                    }
                    else
                    {
                        emissionProb[i][t] = 0;
                    }
                }
                // no flags, no strand-specific observations
                else                              
                {
                    emissionProb[i][t] = this->emissions[i]->calcEmissionProbability(obs[n][t], isNaN[n][t]);
                }

            }

        }
    }

    free(rev_obs);
}


void HMM::getAlphaBeta(double*** obs, double** alpha, double** beta, double* c, double** emissionProb, int* T, int n, int ncores, double effective_zero, int verbose)
{

    int i,j,k,t;
    int nzero_all = 0;                            //0;
    double **myTransMat = this->A->getTransMat();

    int* nnonzeros = (int*)malloc(sizeof(int)*K);
    int** transitions = (int**)malloc(sizeof(int*)*K);
    int** transitions_rev = (int**)malloc(sizeof(int*)*K);
    int* nnonzeros_rev = (int*)malloc(sizeof(int)*K);
    for(i=0; i<K; i++)
    {
        int currnonzeros = 0;
        int currnonzeros_rev = 0;
        for(j=0; j<K; j++)
        {
            if(this->A->getTransMat()[i][j] > effective_zero)
            {
                currnonzeros++;
            }
            else
            {
                nzero_all++;
            }

            if(this->A->getTransMat()[j][i] > effective_zero)
            {
                currnonzeros_rev++;
            }
        }
        nnonzeros[i] = currnonzeros;
        transitions[i] = (int*)malloc(sizeof(int)*currnonzeros);
        nnonzeros_rev[i] = currnonzeros_rev;
        transitions_rev[i] = (int*)malloc(sizeof(int)*currnonzeros_rev);
        currnonzeros = 0;
        currnonzeros_rev = 0;
        for(j=0; j<K; j++)
        {
            if(myTransMat[i][j] > effective_zero)
            {
                transitions[i][currnonzeros++] = j;
            }
            if(myTransMat[j][i] > effective_zero)
            {
                transitions_rev[i][currnonzeros_rev++] = j;
            }
        }
    }

    if(verbose)
    {
        Rprintf("Sequence %d => calculating forward and backward terms (%d transitions are effectively 0).                                      \r", n+1, nzero_all);
    }


    // Compute alpha_1(i) and rescale with c_1
    c[0] = 0;
    for(i=0; i<this->K; i++)
    {
        alpha[0][i] = this->pi->getInitialProb()[i]*emissionProb[i][0];
        c[0] = c[0] + alpha[0][i];
    }
    c[0] = 1/c[0];
    for(i=0; i<this->K; i++)
    {
        alpha[0][i] = c[0]*alpha[0][i];
        if(alpha[0][i] < 1e-300)
        {
            alpha[0][i] = 1e-300;
        }
    }

    // Compute all other alpha_t(i)
    for(t=1; t<T[n]; t++)
    {
        R_CheckUserInterrupt();
        c[t] = 0;
        for(i=0; i<K; i++)
        {
            alpha[t][i] = 0;

            for(k=0; k<nnonzeros_rev[i]; k++)
            {
                alpha[t][i] = alpha[t][i] + alpha[t-1][transitions_rev[i][k]]*myTransMat[transitions_rev[i][k]][i];
            }
            alpha[t][i] = alpha[t][i]*emissionProb[i][t];
            if(alpha[t][i] < 1e-300)
            {
                alpha[t][i] = 1e-300;             //effective_zero;
            }
            c[t] = c[t] + alpha[t][i];

        }
        // Rescale alpha[t][i]
        c[t] = 1/c[t];

        for(i=0; i<K; i++)
        {
            alpha[t][i] = c[t]*alpha[t][i];
        }
    }

    // beta_T = 1 is rescaled
    for(i=0; i<K; i++)
    {
        beta[T[n]-1][i] = 1;
    }

    // Compute all other beta_t(i)
    for(t=T[n]-2; t>=0; t--)
    {
        R_CheckUserInterrupt();
        for(i=0; i<K; i++)
        {
            beta[t][i] = 0;
            for(k=0; k<nnonzeros[i]; k++)
            {
                beta[t][i] = beta[t][i]+myTransMat[i][transitions[i][k]]*emissionProb[transitions[i][k]][t+1]*beta[t+1][transitions[i][k]];
            }
            if(beta[t][i] < 1e-300)
            {
                beta[t][i] = 1e-300;
            }
            beta[t][i] = c[t]*beta[t][i];
        }
    }


}


void HMM::getGamma(double** alpha, double** beta, double* c,  double** emissionProb, double** gamma, int* T, int n, int ncores, double effective_zero, int verbose)
{
    int iter;

    int p,q;
    int nzero_all = 0;
    int* nnonzeros = (int*)malloc(sizeof(int)*K);
    int** transitions = (int**)malloc(sizeof(int*)*K);
    for(p=0; p<K; p++)
    {
        int currnonzeros = 0;
        for(q=0; q<K; q++)
        {
            if(this->A->getTransMat()[p][q] > effective_zero)
            {
                currnonzeros++;
            }
            else
            {
                nzero_all++;
            }
        }
        nnonzeros[p] = currnonzeros;
        transitions[p] = (int*)malloc(sizeof(int)*currnonzeros);
        currnonzeros = 0;
        for(q=0; q<K; q++)
        {
            if(this->A->getTransMat()[p][q] > effective_zero)
            {
                transitions[p][currnonzeros++] = q;
            }
        }
    }

    if(verbose)
    {
        Rprintf("Sequence %d => calculating gamma (%d transitions are effectively 0).                                             \r", n+1, nzero_all);
    }

    int seglen = (int) T[n]/ncores;
    int* myBounds = (int*)malloc(sizeof(int)*(ncores+1));

    myBounds[0] = 0;
    for(iter=1; iter<ncores; iter++)
    {
        myBounds[iter] = myBounds[iter-1]+seglen;
    }
    myBounds[ncores] = T[n];

    int curr_core;
    double **myTransMat = this->A->getTransMat();

    #pragma omp parallel for num_threads(ncores)
    for(curr_core=0; curr_core<ncores; curr_core++)
    {
        int i,j,t,k;
        double denom;
        for(t=myBounds[curr_core]; t<myBounds[curr_core+1]; t++)
        {

            // gamma
            denom = 0.0;
            for(i=0 ; i<this->K; i++)
            {
                gamma[t][i] = alpha[t][i] * beta[t][i];
                denom += gamma[t][i];
            }
            for (i=0; i<this->K ; i++)
            {
                gamma[t][i] /= denom ;
            }
        }
    }

}


void HMM::getGammaXsi(double** alpha, double** beta, double* c,  double** emissionProb, double** gamma, double*** xsi, int* T, int n, int ncores, double effective_zero, int verbose)
{
    int iter;

    int p,q;
    int nzero_all = 0;
    int* nnonzeros = (int*)malloc(sizeof(int)*K);
    int** transitions = (int**)malloc(sizeof(int*)*K);
    for(p=0; p<K; p++)
    {
        int currnonzeros = 0;
        for(q=0; q<K; q++)
        {
            if(this->A->getTransMat()[p][q] > effective_zero)
            {
                currnonzeros++;
            }
            else
            {
                nzero_all++;
            }
        }
        nnonzeros[p] = currnonzeros;
        transitions[p] = (int*)malloc(sizeof(int)*currnonzeros);
        currnonzeros = 0;
        for(q=0; q<K; q++)
        {
            if(this->A->getTransMat()[p][q] > effective_zero)
            {
                transitions[p][currnonzeros++] = q;
            }
        }
    }

    if(verbose)
    {
        Rprintf("Sequence %d => calculating gamma and xi terms  (%d transitions are effectively 0).                                             \r", n+1, nzero_all);
    }

    int seglen = (int) T[n]/ncores;
    int* myBounds = (int*)malloc(sizeof(int)*(ncores+1));

    myBounds[0] = 0;
    for(iter=1; iter<ncores; iter++)
    {
        myBounds[iter] = myBounds[iter-1]+seglen;
    }
    myBounds[ncores] = T[n];

    int curr_core;
    double **myTransMat = this->A->getTransMat();

    #pragma omp parallel for num_threads(ncores)
    for(curr_core=0; curr_core<ncores; curr_core++)
    {
        int i,j,t,k;
        double denom;
        for(t=myBounds[curr_core]; t<myBounds[curr_core+1]; t++)
        {

            // gamma
            denom = 0.0;
            for(i=0 ; i<this->K; i++)
            {
                gamma[t][i] = alpha[t][i] * beta[t][i];
                denom += gamma[t][i];
            }
            for (i=0; i<this->K ; i++)
            {
                gamma[t][i] /= denom ;
            }

            // xsi
            if(t<T[n]-1)
            {
                for(i=0 ; i<this->K; i++)
                {
                    denom = 1/c[t] * beta[t][i];
                    /*for(k=0 ; k<this->K; k++) {
                        xsi[t][i][k] = (gamma[t][i]*myTransMat[i][k]*emissionProb[k][t+1] * beta[t+1][k])  / denom;
                    }*/
                    for(k=0 ; k<nnonzeros[i]; k++)
                    {
                        xsi[t][i][transitions[i][k]] = (gamma[t][i]*myTransMat[i][transitions[i][k]]*emissionProb[transitions[i][k]][t+1] * beta[t+1][transitions[i][k]])  / denom;
                    }
                }
            }
        }
    }

}


int HMM::allocateMemEM(double*** emissionProb, double*** alpha, double*** beta, double*** gamma, double**** xsi, double** c, double** Pk, int maxLen, int nsample)
{
    int i,j,t;

    int memory_used = 0;

    *c = (double*)malloc(sizeof(double)*maxLen);
    memory_used += sizeof(double)*maxLen;
    *emissionProb = (double**)malloc(sizeof(double*)*this->K);
    memory_used += sizeof(double*)*this->K;

    for(i=0; i<this->K; i++)
    {
        (*emissionProb)[i] = (double*)malloc(sizeof(double)*maxLen);
        memory_used += sizeof(double)*maxLen;
        for(t=0; t<maxLen; t++)
        {
            (*emissionProb)[i][t] = 0;
        }
    }

    *alpha = (double**)malloc(sizeof(double*)*maxLen);
    memory_used += sizeof(double*)*maxLen;
    *beta  = (double**)malloc(sizeof(double*)*maxLen);
    memory_used += sizeof(double*)*maxLen;
    *gamma  = (double**)malloc(sizeof(double*)*maxLen);
    memory_used += sizeof(double*)*maxLen;
    *xsi  = (double***)malloc(sizeof(double**)*maxLen);
    memory_used += sizeof(double**)*maxLen;

    for(t=0; t<maxLen; t++)
    {
        (*c)[t] = 0;
        (*alpha)[t] = (double*)malloc(sizeof(double)*this->K);
        memory_used += sizeof(double)*this->K;
        (*beta)[t] = (double*)malloc(sizeof(double)*this->K);
        memory_used += sizeof(double)*this->K;
        (*gamma)[t] = (double*)malloc(sizeof(double)*this->K);
        memory_used += sizeof(double)*this->K;
        (*xsi)[t] = (double**)malloc(sizeof(double*)*this->K);
        memory_used += sizeof(double*)*this->K;

        for(i=0; i<K; i++)
        {
            (*alpha)[t][i] = 0;
            (*beta)[t][i] = 0;
            (*gamma)[t][i] = 0;
            (*xsi)[t][i] = (double*)malloc(sizeof(double)*this->K);
            memory_used += sizeof(double)*this->K;
            for(j=0; j<K; j++)
            {
                (*xsi)[t][i][j] = 0;
            }
        }
    }

    *Pk = (double*)malloc(sizeof(double)*nsample);
    memory_used += sizeof(double)*nsample;
    double megabytes_used = ((double)memory_used)/1000000;

    if(DEBUG_MEMORY)
    {
        Rprintf("Baum-Welch needs %lf MB of memory.\n", megabytes_used);
    }
    return memory_used;
}


int HMM::deallocateMemEM(double** emissionProb, double** alpha, double** beta, double** gamma, double*** xsi, double* c, double* Pk, int maxLen, int nsample)
{
    int i,j,t;
    int memory_free = 0;

    for(i=0; i<this->K; i++)
    {
        free(emissionProb[i]);
        memory_free += sizeof(double)*maxLen;
    }
    free(emissionProb);
    memory_free += sizeof(double*)*this->K;
    free(c);
    memory_free += sizeof(double)*maxLen;
    free(Pk);
    memory_free += sizeof(double)*nsample;

    for(t=0; t<maxLen; t++)
    {
        free(alpha[t]);
        memory_free += sizeof(double)*this->K;
        free(beta[t]);
        memory_free += sizeof(double)*this->K;
        free(gamma[t]);
        memory_free += sizeof(double)*this->K;

        for(i=0; i<K; i++)
        {
            free(xsi[t][i]);
            memory_free += sizeof(double)*this->K;
        }
        free(xsi[t]);
        memory_free += sizeof(double)*this->K;
    }

    free(alpha);
    memory_free += sizeof(double*)*maxLen;
    free(beta);
    memory_free += sizeof(double*)*maxLen;
    free(gamma);
    memory_free += sizeof(double*)*maxLen;
    free(xsi);
    memory_free += sizeof(double**)*maxLen;

    return memory_free;
}


void HMM::updateSampleAux(double*** observations, int* T, int n, double** alpha, double** beta, double** gamma, double*** xsi, double* Pk, int* state2flag, int* couples, int* revop, int** isNaN, double*** fixedEmission, SEXP bidirOptimParams, SEXP emissionPrior, int ncores, double effective_zero, int verbose)
{
    int i;
    // udpate auxiliary terms for initial state probabilities
    for(i=0; i<this->K; i++)
    {
        if(couples == NULL)
        {
            this->pi->updateSample(gamma, i);
        }
        else
        {
            this->pi->updateSampleCoupled(gamma, i, couples, bidirOptimParams);
        }
    }

    // update auxiliary terms for transition matrix
    if(couples == NULL)
    {
        this->A->updateAuxiliaries(gamma, xsi,  Pk, T, n, isNaN, ncores, effective_zero, verbose);

    }
    else
    {
        this->A->updateAuxiliaries(gamma, xsi,  Pk, T, n, couples, bidirOptimParams, isNaN, ncores, effective_zero, verbose);
    }

    if(this->K < ncores)
    {
        ncores = this->K;
    }
    int *myStateBuckets = (int*)malloc(sizeof(int)*ncores+1);
    for(i=0; i<=ncores; i++)
    {
        myStateBuckets[i] = 0;
    }
    int currbucket = 1;
    for(i=0; i<this->K; i++)
    {
        myStateBuckets[currbucket] = myStateBuckets[currbucket]+1;
        if(currbucket == ncores)
        {
            currbucket = 0;
        }
        currbucket++;
    }
    for(i=1; i<ncores+1; i++)
    {
        myStateBuckets[i] = myStateBuckets[i]+myStateBuckets[i-1];
    }

    // update auxiliary terms for emission functions
    if(fixedEmission == NULL)
    {

        if(verbose)
        {
            Rprintf("Sequence %d => Updating emission distributions auxiliary terms.                                     \r", n+1);
        }

        int k;
        #pragma omp parallel for num_threads(ncores)
        for(k=1; k<ncores+1; k++)
        {

            int state;
            for(state=myStateBuckets[k-1]; state<myStateBuckets[k]; state++)
            {
                // coupled states and reverse observations
                if(couples != NULL && revop != NULL && couples[state] != state)
                {
                    if(DEBUG)
                    {
                        Rprintf("Updating emission auxiliaries: +couples +rev.observations\n");
                    }
                    int p = couples[state];
                    this->emissions[state]->updateAuxiliariesCoupledRevop(observations, gamma, Pk, T, n, state, p, state2flag, revop, isNaN);
                }
                // coupled states and no reverse observations
                else if(couples != NULL && revop == NULL && state2flag != NULL && couples[state] != state)
                {
                    if(DEBUG)
                    {
                        Rprintf("Updating emission auxiliaries: +couples -rev.observations\n");
                    }
                    int p = couples[state];
                    this->emissions[state]->updateAuxiliariesCoupled(observations, gamma, Pk, T, n, state, p, isNaN);
                }
                // no coupled states and no reverse observations
                else                              
                {
                    if(DEBUG)
                    {
                        Rprintf("Updating emission auxiliaries: -couples -rev.observations\n");
                    }
                    this->emissions[state]->updateAuxiliaries(observations, gamma, Pk, T, n, state, isNaN);
                }
            }
        }
    }
}


void HMM::updateModelParams(double*** observations, int nsample, int* state2flag, int* couples, int* revop, int verbose, int updateTransMat, int** isNaN, double*** fixedEmission, SEXP bidirOptimParams, SEXP emissionPrior, int ncores, double effective_zero, int* myStateBuckets, double* Pk, int curriter, int currN)
{

    // Update transition matrix
    if(updateTransMat == 1)
    {
        if(LENGTH(bidirOptimParams) > 0)
        {
            this->A->update(bidirOptimParams);
        }
        else if(couples == NULL && updateTransMat == 1)
        {
            this->A->update(effective_zero);
        }
        else
        {
            this->A->update(couples, effective_zero);
        }
    }

    // update initial state probabilities
    this->pi->update(nsample, bidirOptimParams);

    // update parameters of emissions
    if(fixedEmission == NULL)
    {
        if(verbose)
        {
            Rprintf("Updating emission distributions.                                                            \r");
        }
        int k;
        #pragma omp parallel for num_threads(ncores)
        for(k=1; k<ncores+1; k++)
        {
            int state;
            for(state=myStateBuckets[k-1]; state<myStateBuckets[k]; state++)
            {
                // coupled states and reverse observations
                if(couples != NULL && revop != NULL && couples[state] != state)
                {
                    int p = couples[state];
                    if(DEBUG)
                    {
                        Rprintf("Updating emission parameters: +couples +rev.observations\n");
                    }
                    this->emissions[state]->updateCoupledRevop(observations, Pk, p, state2flag, revop, this->emissions[p]->getParameter()->getGammaAux(), isNaN, emissionPrior, currN);
                }
                // coupled states and no reverse observations
                else if(couples != NULL && revop == NULL && state2flag != NULL && couples[state] != state)
                {
                    if(DEBUG)
                    {
                        Rprintf("Updating emission parameters: +couples -rev.observations\n");
                    }
                    this->emissions[state]->update(observations, Pk, isNaN, emissionPrior, currN);
                }
                else                              // no coupled states and no reverse observations
                {
                    if(DEBUG)
                    {
                        Rprintf("Updating emission parameters: -couples -rev.observations\n");
                    }
                    this->emissions[state]->update(observations, Pk, isNaN, emissionPrior, currN);
                }
            }
        }
    }
}


list<double> HMM::BaumWelch(double*** observations, int* T, int nsample, int maxIters, int** flags, int* state2flag, int* couples, int* revop, int verbose, int updateTransMat, int** isNaN, double*** fixedEmission, SEXP bidirOptimParams, SEXP emissionPrior, int ncores, double effective_zero, double convergence, int incrementalEM)
{
    int t,i,j,n;

    list<double> log_lik;
    double old_prior =  -(double)INFINITY;
    int iter = 0;
    double old_log_lik = -(double)INFINITY;
    double new_log_lik = 0;

    int maxLen = 0;
    for(n=0; n<nsample; n++)
    {
        if(maxLen < T[n])
        {
            maxLen = T[n];
        }
    }

    if(this->K < ncores)
    {
        ncores = this->K;
    }
    int *myStateBuckets = (int*)malloc(sizeof(int)*ncores+1);
    for(i=0; i<=ncores; i++)
    {
        myStateBuckets[i] = 0;
    }
    int currbucket = 1;
    for(i=0; i<this->K; i++)
    {
        myStateBuckets[currbucket] = myStateBuckets[currbucket]+1;
        if(currbucket == ncores)
        {
            currbucket = 0;
        }
        currbucket++;
    }
    for(i=1; i<ncores+1; i++)
    {
        myStateBuckets[i] = myStateBuckets[i]+myStateBuckets[i-1];
    }

    // allocate memory for auxiliary variables
    double*** xsi = NULL;
    double** alpha = NULL;
    double** beta = NULL;
    double** gamma = NULL;
    double* Pk = NULL;
    double** emissionProb = NULL;
    double* c = NULL;
    int memory_used = this->allocateMemEM(&emissionProb, &alpha, &beta, &gamma, &xsi, &c, &Pk, maxLen, nsample);
    const char* constraints_changed = "";

    int before = 0;
    int now = 0;
    double currDiff = (double)INFINITY;

    // calculate Log-Likelihood and E-Step for inital model
    for(n=0; n<nsample; n++)
    {

        // calculate emission probabilities
        if(fixedEmission == NULL)
        {
            this->calcEmissionProbs(observations, emissionProb, T, n, flags, state2flag, revop, isNaN, ncores, verbose);
        }
        else
        {
            if(fixedEmission != NULL)
            {
                for(i=0; i<this->K; i++)
                {
                    for(t=0; t<T[n]; t++)
                    {
                        emissionProb[i][t] = fixedEmission[n][i][t];
                    }
                }
            }
        }
        // calculate alpha and beta
        this->getAlphaBeta(observations, alpha, beta, c, emissionProb, T, n, ncores, effective_zero, verbose);

        // Calculate gamma and xsi
        this->getGammaXsi(alpha, beta, c, emissionProb, gamma, xsi, T, n, ncores, effective_zero, verbose);

        // Calculate weight of the current sample
        Pk[n] = 0;
        for(i=0; i<this->K; i++)
        {
            Pk[n] += alpha[T[n]-1][i];
        }

        // calculate Log-Likelihood of the current model
        for(t=0; t<T[n]; t++)
        {
            if(isNaN[n][t] == 0)
            {
                new_log_lik = new_log_lik + log(c[t]);
            }
        }

        // update auxiliary terms for the current sample n
        this->updateSampleAux(observations, T, n, alpha, beta, gamma, xsi, Pk, state2flag, couples, revop, isNaN, fixedEmission, bidirOptimParams, emissionPrior, ncores, effective_zero, verbose);
        if(incrementalEM)
        {
            this->updateModelParams(observations, nsample, state2flag, couples, revop, verbose, updateTransMat, isNaN, fixedEmission, bidirOptimParams, emissionPrior, ncores, effective_zero, myStateBuckets, Pk, iter, n);
        }
    }

    // add prior for emission (gaussian)
    const char* llh = "Log-Likelihood";
    double prior = 0;

    if(LENGTH(emissionPrior) != 0)
    {
        for(i=0; i<K; i++)
        {
            double temp_prior = this->emissions[i]->Prior(emissionPrior);
            prior += temp_prior;
        }
        old_prior = prior;
        llh = "Log-Posterior";
    }

    new_log_lik = -new_log_lik+prior;
    if(verbose)
    {
        Rprintf("Iteration %d: %s of initial model: %f                                                        \n", 0, llh, new_log_lik);
    }

    while(iter == 0 | ((iter < maxIters &&  fabs(currDiff)>convergence) | (iter < maxIters &&  now != before) ))
    {
        R_CheckUserInterrupt();
        old_log_lik = new_log_lik;
        new_log_lik = 0;

        if(! incrementalEM)
        {
            this->updateModelParams(observations, nsample, state2flag, couples, revop, verbose, updateTransMat, isNaN, fixedEmission, bidirOptimParams, emissionPrior, ncores, effective_zero, myStateBuckets, Pk, iter, -1);
        }

        // begin new sample (e.g. chromosome) => calculate E-Step for all samples of new model
        for(n=0; n<nsample; n++)
        {

            // calculate emission probabilities
            if(fixedEmission == NULL)
            {
                this->calcEmissionProbs(observations, emissionProb, T, n, flags, state2flag, revop, isNaN, ncores, verbose);
            }
            else
            {
                if(fixedEmission != NULL)
                {
                    for(i=0; i<this->K; i++)
                    {
                        for(t=0; t<T[n]; t++)
                        {
                            emissionProb[i][t] = fixedEmission[n][i][t];
                        }
                    }
                }
            }

            // calculate alpha and beta
            this->getAlphaBeta(observations, alpha, beta, c, emissionProb, T, n, ncores, effective_zero, verbose);

            // Calculate gamma and xsi
            this->getGammaXsi(alpha, beta, c, emissionProb, gamma, xsi, T, n, ncores, effective_zero, verbose);

            // Calculate weight of the current sample
            Pk[n] = 0;
            for(i=0; i<this->K; i++)
            {
                Pk[n] += alpha[T[n]-1][i];
            }

            // calculate Log-Likelihood of the current model
            for(t=0; t<T[n]; t++)
            {
                if(isNaN[n][t] == 0)
                {
                    new_log_lik = new_log_lik + log(c[t]);
                }
            }

            // calculate new LLH
            if(LENGTH(bidirOptimParams) > 0 & iter > 1)
            {
                before = INTEGER(getListElement(bidirOptimParams, "nrm"))[LENGTH(getListElement(bidirOptimParams, "nrm"))-3];
                now = INTEGER(getListElement(bidirOptimParams, "nrm"))[LENGTH(getListElement(bidirOptimParams, "nrm"))-2];
            }

            // update auxiliary terms for the current sample n
            this->updateSampleAux(observations, T, n, alpha, beta, gamma, xsi, Pk, state2flag, couples, revop, isNaN, fixedEmission, bidirOptimParams, emissionPrior, ncores, effective_zero, verbose);
            if(incrementalEM)
            {
                this->updateModelParams(observations, nsample, state2flag, couples, revop, verbose, updateTransMat, isNaN, fixedEmission, bidirOptimParams, emissionPrior, ncores, effective_zero, myStateBuckets, Pk, iter, n);
            }
        }
        // end current sample

        // add prior for emission (gaussian)
        const char* llh = "Log-Likelihood";
        double prior = 0;
        if(LENGTH(emissionPrior) != 0)
        {
            for(i=0; i<K; i++)
            {
                double temp_prior = this->emissions[i]->Prior(emissionPrior);
                prior += temp_prior;
            }
            old_prior = prior;
            llh = "Log-Posterior";
        }

        new_log_lik = -new_log_lik+prior;
        currDiff = (new_log_lik-old_log_lik)/fabs(old_log_lik);
        constraints_changed = "";
        if(iter > 0 && LENGTH(bidirOptimParams) != 0 & now  != before & DEBUG)
        {
            constraints_changed = " (*)";
        }

        iter++;
        if(verbose)
        {
            Rprintf("Iteration %d: %s=%f (relative Diff=%f)%s                                             \n", iter, llh, new_log_lik, fabs(currDiff), constraints_changed);
        }
        log_lik.push_back(new_log_lik);

    }
    // end while

    int memory_free = this->deallocateMemEM(emissionProb, alpha, beta, gamma, xsi, c, Pk, maxLen, nsample);
    if(memory_free == memory_used && DEBUG_MEMORY)
    {
        Rprintf("Memory for auxiliary Baum-Welch variables successfully deallocated.\n");
    }

    return log_lik;
}


void HMM::Viterbi(int **S, double*** obs, int nsample, int* T, int verbose, int** isNaN, double*** fixedEmission)
{

    if(verbose)
    {
        Rprintf("Calculating Viterbi path.\n");
    }

    int i,j,t,n;

    for(n=0; n<nsample; n++)
    {
        R_CheckUserInterrupt();
        double **delta = (double**)malloc(sizeof(double*)*T[n]);
        int **psi = (int**)malloc(sizeof(int*)*T[n]);

        for(t=0; t<T[n]; t++)
        {
            delta[t] = (double*)malloc(sizeof(double)*this->K);
            psi[t] = (int*)malloc(sizeof(int)*this->K);
        }
        double curr_emission;
        // Initialization
        for(i=0; i<this->K; i++)
        {
            if(fixedEmission == NULL)
            {
                delta[0][i] = log(this->pi->getInitialProb()[i])+log(this->emissions[i]->calcEmissionProbability(obs[n][0], isNaN[n][0]));
            }
            else
            {
                curr_emission = fixedEmission[n][i][0];
                if(curr_emission < 1e-100)
                {
                    curr_emission = 1e-100;
                }

                delta[0][i] = log(this->pi->getInitialProb()[i])+log(curr_emission);
            }
            psi[0][i] = 0;
        }

        // Recursion
        for(t=1; t<T[n]; t++)
        {
        R_CheckUserInterrupt();
            for(j=0; j<this->K; j++)
            {
                delta[t][j] = -(double)INFINITY;
                int argmaxi = -1;
                double curr_max = -(double)INFINITY;
                for(i=0; i<this->K; i++)
                {
                    double curr_val = 0;

                    if(fixedEmission == NULL)
                    {
                        curr_val = delta[t-1][i]+log(this->A->getTransMat()[i][j])+log(this->emissions[j]->calcEmissionProbability(obs[n][t], isNaN[n][t]));
                    }
                    else
                    {
                        curr_emission = fixedEmission[n][j][t];
                        if(curr_emission < 1e-100)
                        {
                            curr_emission = 1e-100;
                        }
                        curr_val = delta[t-1][i]+log(this->A->getTransMat()[i][j])+log(curr_emission);
                    }
                    if(curr_val > delta[t][j])
                    {
                        delta[t][j] = curr_val;
                    }
                    if(curr_val > curr_max)
                    {
                        curr_max = curr_val;
                        argmaxi = i;
                    }
                }
                psi[t][j] = argmaxi;
            }
        }

        // Termination
        double P = -(double)INFINITY;
        for(i=0; i<this->K; i++)
        {
            if(P < delta[T[n]-1][i])
            {
                P = delta[T[n]-1][i];
                S[n][T[n]-1] = i;
            }
        }

        for(t=T[n]-2; t>=0; t--)
        {
            S[n][t] = psi[t+1][S[n][t+1]];
        }

        for(t=0; t<T[n]; t++)
        {
            free(delta[t]);
            free(psi[t]);
        }

        free(delta);
        free(psi);
        if(verbose)
        {
            Rprintf("Viterbi path #%d. LLH=%f\n", n+1, P);
        }
    }
}
