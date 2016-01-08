
#include "TransitionMatrix.h"

TransitionMatrix::TransitionMatrix(double **A, int K)
{
    this->A = A;
    this->K = K;

    this->updateNumerator = (double**)malloc(sizeof(double*)*K);
    this->updateDenominator = (double**)malloc(sizeof(double*)*K);
    int i,j;
    for(i=0; i<K; i++)
    {
//printf("i=%d\n", i);
        this->updateNumerator[i] = (double*)malloc(sizeof(double)*K);
        this->updateDenominator[i] = (double*)malloc(sizeof(double)*K);
        for(j=0; j<K; j++)
        {
//	printf("1.1\n");
            this->updateNumerator[i][j] = 0;
            this->updateDenominator[i][j] = 0;
//	printf("1.2\n");
        }
    }
    int mem = 3*sizeof(double*)*this->K + 3*sizeof(double)*this->K;
    if(DEBUG_MEMORY)
    {
        printf("new->TransitionMatrix; (%d bytes)\n", mem);
    }
}


TransitionMatrix::~TransitionMatrix()
{
    int i;
    for(i=0; i<this->K; i++)
    {
        free(this->A[i]);
        free(this->updateNumerator[i]);
        free(this->updateDenominator[i]);
    }
    free(this->A);
    free(this->updateNumerator);
    free(this->updateDenominator);

    int mem = 3*sizeof(double*)*this->K + 3*sizeof(double)*this->K;
    if(DEBUG_MEMORY)
    {
        Rprintf("delete->TransitionMatrix; (%d bytes)\n", mem);
    }
}


double** TransitionMatrix::getTransMat()
{
    return this->A;
}


void TransitionMatrix::updateAuxiliaries(double** gamma, double*** xsi,  double* Pk, int* T, int n, int* couples, SEXP bidirOptimParams, int** isNaN, int ncores, double effective_zero, int verbose)
{
//effective_zero=-1;
    int a,b;
    int nzero_all = 0;
    int* nnonzeros = (int*)malloc(sizeof(int)*K);
    int** transitions = (int**)malloc(sizeof(int*)*K);
    for(a=0; a<K; a++)
    {
        int currnonzeros = 0;
        for(b=0; b<K; b++)
        {
            if(this->A[a][b] > effective_zero)
            {
                currnonzeros++;
            }
            else
            {
                nzero_all++;
            }
        }
        nnonzeros[a] = currnonzeros;
        transitions[a] = (int*)malloc(sizeof(int)*currnonzeros);
        currnonzeros = 0;
        for(b=0; b<K; b++)
        {
            if(this->A[a][b] > effective_zero)
            {
                transitions[a][currnonzeros++] = b;
            }
        }
    }

    int i,j,k;
    double numer, denom;
    double* xsi_sum = (double*)malloc(sizeof(double)*this->K);
    if(verbose)
    {
        Rprintf("Sequence %d => Updating transition auxiliaries (%d transitions are effectively 0).                                     \r", n+1, nzero_all);
    }

    if(LENGTH(bidirOptimParams) > 0)
    {
        int t;
        for(i=0; i<this->K; i++)
        {
            xsi_sum[i] = 0.0;
            for(j=0; j<K; j++)
            {
                numer = 0.0;
                denom = 0.0;
                for(t=1; t<T[n]; t++)
                {
//if(isNaN[n][t] == 0 & isNaN[n][t-1] == 0) {
                    numer = numer + xsi[t-1][i][j];
//}
                }
// 				// the numerator holds sum(xsi) over all samples for numerical optimization using Rsolnp
                                                  //1/Pk[n]*
                this->updateNumerator[i][j] += numer;
                this->updateDenominator[i][j] = 0;
            }
        }
    }
    else                                          // calculate analytical update for bdHMM transitions
    {

//parallelizes over ncores
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
/*for(i=0; i<=ncores; i++) {
    printf("bucket[%d]=%d\n", i, myStateBuckets[i]);
}*/
// parallelization over states done

        int *twin_of = (int*)malloc(sizeof(int)*this->K);
        for(i=0; i<this->K; i++)
        {
            twin_of[i] = couples[i];
        }

        int i,j;
        double** numer = (double**)malloc(sizeof(double*)*this->K);
        double* denom = (double*)malloc(sizeof(double)*this->K);
        for(i=0; i<this->K; i++)
        {
            numer[i] = (double*)malloc(sizeof(double)*this->K);
            denom[i] = 0;
        }

        for(i=0; i<this->K; i++)
        {
            for(j=0; j<K; j++)
            {
                numer[i][j] = 0.0;
            }
        }

        int k;
#pragma omp parallel for
        for(k=1; k<ncores+1; k++)
        {
            int p;
//printf("bucket[%d]=(%d:%d)\n", k, myStateBuckets[k-1], myStateBuckets[k]-1);
            for(p=myStateBuckets[k-1]; p<myStateBuckets[k]; p++)
            {
                int q,t;
                for(t=1; t<T[n]; t++)
                {
                    for(q=0 ; q<nnonzeros[p]; q++)
                    {
                        numer[p][transitions[p][q]] = numer[p][transitions[p][q]] + (xsi[t-1][twin_of[transitions[p][q]]][twin_of[p]] + xsi[t-1][p][transitions[p][q]]);
                    }
                    denom[p] = denom[p] + (gamma[t][twin_of[p]] + gamma[t-1][p]);

                }
            }
        }

        for(i=0; i<this->K; i++)
        {
            for(j=0; j<K; j++)
            {
                this->updateNumerator[i][j] += 1/Pk[n]*numer[i][j];
                this->updateDenominator[i][j] += 1/Pk[n]*denom[i];
            }
        }

        for(i=0; i<this->K; i++)
        {
            free(numer[i]);
        }
        free(numer);
        free(denom);
        free(myStateBuckets);
        free(twin_of);
    }

    free(nnonzeros);
    for(i=0; i<K; i++)
    {
        free(transitions[i]);
    }
    free(transitions);

    free(xsi_sum);
}


void TransitionMatrix::updateAuxiliaries(double** gamma, double*** xsi,  double* Pk, int* T, int n, int** isNaN, int ncores, double effective_zero, int verbose)
{

    int a,b;
    int nzero_all = 0;
    int* nnonzeros = (int*)malloc(sizeof(int)*K);
    int** transitions = (int**)malloc(sizeof(int*)*K);
    for(a=0; a<K; a++)
    {
        int currnonzeros = 0;
        for(b=0; b<K; b++)
        {
            if(this->A[a][b] > effective_zero)
            {
                currnonzeros++;
            }
            else
            {
                nzero_all++;
            }
        }
        nnonzeros[a] = currnonzeros;
        transitions[a] = (int*)malloc(sizeof(int)*currnonzeros);
        currnonzeros = 0;
        for(b=0; b<K; b++)
        {
            if(this->A[a][b] > effective_zero)
            {
                transitions[a][currnonzeros++] = b;
            }
        }
    }

    int i,j;
//parallelizes over ncores
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
/*for(i=0; i<=ncores; i++) {
    printf("bucket[%d]=%d\n", i, myStateBuckets[i]);
}*/
// parallelization over states done

    double** numer = (double**)malloc(sizeof(double*)*this->K);
    double* denom = (double*)malloc(sizeof(double)*this->K);
    for(i=0; i<this->K; i++)
    {
        numer[i] = (double*)malloc(sizeof(double)*this->K);
        denom[i] = 0;
    }

    for(i=0; i<this->K; i++)
    {
        for(j=0; j<K; j++)
        {
            numer[i][j] = 0.0;
        }
    }

    if(verbose)
    {
        Rprintf("Sequence %d => Updating transition auxiliaries.                                     \r", n+1);
    }

    int k;
#pragma omp parallel for
    for(k=1; k<ncores+1; k++)
    {
        int p;
//	printf("bucket[%d]=(%d:%d)\n", k, myStateBuckets[k-1], myStateBuckets[k]-1);
        for(p=myStateBuckets[k-1]; p<myStateBuckets[k]; p++)
        {
            int q,t;
            for(t=1; t<T[n]; t++)
            {
                for(q=0 ; q<nnonzeros[p]; q++)
                {
                    numer[p][transitions[p][q]] = numer[p][transitions[p][q]] + xsi[t-1][p][transitions[p][q]];
                }
                denom[p] = denom[p] + gamma[t-1][p];
            }
        }
    }

    for(i=0; i<this->K; i++)
    {
        for(j=0; j<K; j++)
        {
            this->updateNumerator[i][j] += 1/Pk[n]*numer[i][j];
            this->updateDenominator[i][j] += 1/Pk[n]*denom[i];
        }
    }

    for(i=0; i<this->K; i++)
    {
        free(numer[i]);
    }
    free(numer);
    free(denom);

    free(nnonzeros);
    for(i=0; i<K; i++)
    {
        free(transitions[i]);
    }
    free(transitions);
    free(myStateBuckets);
}


void TransitionMatrix::update(double effective_zero)
{
    int a,b;
    int nzero_all = 0;
    int* nnonzeros = (int*)malloc(sizeof(int)*K);
    int** transitions = (int**)malloc(sizeof(int*)*K);
    for(a=0; a<K; a++)
    {
        int currnonzeros = 0;
        for(b=0; b<K; b++)
        {
            if(this->A[a][b] > effective_zero)
            {
                currnonzeros++;
            }
            else
            {
                nzero_all++;
                this->A[a][b] = 0;
            }
        }
        nnonzeros[a] = currnonzeros;
        transitions[a] = (int*)malloc(sizeof(int)*currnonzeros);
        currnonzeros = 0;
        for(b=0; b<K; b++)
        {
            if(this->A[a][b] > effective_zero)
            {
                transitions[a][currnonzeros++] = b;
            }
        }
    }

    int i,j,k;
//	Rprintf("here->update\n");
    for(i=0; i<this->K; i++)
    {
//printf("gamma[%d]=%f, ", i, this->updateDenominator[i][0]);
        for(k=0; k<nnonzeros[i]; k++)
        {
            j = transitions[i][k];
            this->A[i][j] = this->updateNumerator[i][j]/this->updateDenominator[i][j];
//printf("%f/%f ", this->updateNumerator[i][j], this->updateDenominator[i][j]);
            this->updateNumerator[i][j] = 0;
            this->updateDenominator[i][j] = 0;
//printf("%f ", this->A[i][j]);
        }
//	printf("\n");
//printf("\n");
    }

    free(nnonzeros);
    for(i=0; i<K; i++)
    {
        free(transitions[i]);
    }
    free(transitions);

}


void TransitionMatrix::update(int* couples, double effective_zero)
{
// effective_zero=-1;
    int a,b;
    int nzero_all = 0;
    int* nnonzeros = (int*)malloc(sizeof(int)*K);
    int** transitions = (int**)malloc(sizeof(int*)*K);
    for(a=0; a<K; a++)
    {
        int currnonzeros = 0;
        for(b=0; b<K; b++)
        {
            if(this->A[a][b] > effective_zero)
            {
                currnonzeros++;
            }
            else
            {
                nzero_all++;
                this->A[a][b] = 0;
            }
        }
        nnonzeros[a] = currnonzeros;
        transitions[a] = (int*)malloc(sizeof(int)*currnonzeros);
        currnonzeros = 0;
        for(b=0; b<K; b++)
        {
            if(this->A[a][b] > effective_zero)
            {
                transitions[a][currnonzeros++] = b;
            }
        }
    }

    int i,j,k;
//	Rprintf("here->update\n");
    for(i=0; i<this->K; i++)
    {
//printf("gamma[%d]=%f, ", i, this->updateDenominator[i][0]);
        for(k=0; k<nnonzeros[i]; k++)
        {
            j = transitions[i][k];
            this->A[i][j] = this->updateNumerator[i][j]/this->updateDenominator[i][j];
//			printf("%f/%f ", this->updateNumerator[i][j], this->updateDenominator[i][j]);
            this->updateNumerator[i][j] = 0;
            this->updateDenominator[i][j] = 0;
//printf("%f ", this->A[i][j]);
        }
//	printf("\n");
//printf("\n");
    }

    free(nnonzeros);
    for(i=0; i<K; i++)
    {
        free(transitions[i]);
    }
    free(transitions);

}


int getListElementIndex(SEXP list, const char *str)
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    int i;
    for ( i = 0; i < length(list); i++ )
        if ( strcmp(CHAR(STRING_ELT(names, i)), str) == 0 )
    {
        elmt = VECTOR_ELT(list, i);
        break;
    }

    if ( elmt == R_NilValue )
        error("%s missing from list", str);

    return i;
}


void TransitionMatrix::update(SEXP bidirOptimParams)
{
    int i,j;
//printf("\nCalling Rsolnp to optimize transitions...\nInitial parameters: ");

    for(i=0; i<LENGTH(getListElement(bidirOptimParams, "x0")); i++)
    {
//	printf("%f ", REAL(getListElement(bidirOptimParams, "x0"))[i]);
    }
//	printf("\n");

// set current xsi_sum
    SEXP XSISUM = getListElement(bidirOptimParams, "xsi_sum");
    for(i=0; i<this->K; i++)
    {
        for(j=0; j<this->K; j++)
        {
            REAL(XSISUM)[i+this->K*j] = this->updateNumerator[i][j];
//			printf("%d ", REAL(XSISUM)[i+this->K*j] == this->updateNumerator[i][j]);
// 		//	printf("%f ", REAL(XSISUM)[i+this->K*j]);
        }
//		printf("\n");
    }

// call solnp from R for optimization
    SEXP call = PROTECT( lang2( getListElement(bidirOptimParams, "c2optimize"), bidirOptimParams ) ) ;
    SEXP res = PROTECT( eval( call, R_GlobalEnv ) ) ;

    SEXP TRANSMAT = getListElement(res, "transMat");
    SEXP statD = getListElement(res, "statD");
    SEXP newx0 = getListElement(res, "x0");
    SEXP doit = getListElement(res, "doit");
    INTEGER(getListElement(bidirOptimParams, "update"))[0] = INTEGER(doit)[0];
//if(INTEGER(doit)[0] == 1) {

    for(i=0; i<this->K; i++)
    {
        for(j=0; j<this->K; j++)
        {
            this->A[i][j] = REAL(TRANSMAT)[i+this->K*j];
            REAL(getListElement(bidirOptimParams, "transMat"))[i+this->K*j] = this->A[i][j];
//	printf("%f ", this->A[i][j]);
        }
//printf("\n");
    }

    for(i=0; i<LENGTH(statD); i++)
    {
        REAL(getListElement(bidirOptimParams, "statD"))[i] = REAL(statD)[i];
//printf("%f ", REAL(getListElement(bidirOptimParams, "statD"))[i]);
    }
//printf("\n");
//	printf("Updated parameters: ");
    for(i=0; i<LENGTH(newx0); i++)
    {
        REAL(getListElement(bidirOptimParams, "x0"))[i] = REAL(newx0)[i];
//		printf("%f ", REAL(getListElement(bidirOptimParams, "x0"))[i]);
    }
//	printf("\n");

    SEXP objective = getListElement(res, "objective");
    SEXP newobjective;
    PROTECT(newobjective = NEW_NUMERIC(LENGTH(objective)+1));
//printf("%f\n", REAL(objective)[0]);

    for(i=0; i<LENGTH(objective); i++)
    {
        NUMERIC_POINTER(newobjective)[i] = REAL(objective)[i];
//printf("%d\n", i);
//printf("%f ", REAL(newobjective)[i]);
    }
    NUMERIC_POINTER(newobjective)[i] = REAL(objective)[0];
//printf("%f \n", REAL(newobjective)[i]);
    SET_ELEMENT(bidirOptimParams, getListElementIndex(bidirOptimParams, "objective"), newobjective);
    UNPROTECT(1);

    SEXP nrm = getListElement(res, "nrm");
    SEXP newnrm;
    PROTECT(newnrm = NEW_INTEGER(LENGTH((getListElement(bidirOptimParams, "nrm")))+1));
//printf("%f\n", REAL(objective)[0]);

    for(i=0; i<LENGTH(getListElement(bidirOptimParams, "nrm")); i++)
    {
        INTEGER_POINTER(newnrm)[i] = INTEGER(getListElement(bidirOptimParams, "nrm"))[i];
//	printf("%d ", INTEGER(newnrm)[i]);
    }
    INTEGER_POINTER(newnrm)[i] = INTEGER(nrm)[0];
//printf("%d \n", INTEGER(newnrm)[i]);
    SET_ELEMENT(bidirOptimParams, getListElementIndex(bidirOptimParams, "nrm"), newnrm);
    UNPROTECT(1);

//
// 	SEXP nrm = getListElement(res, "nrm");
// 	//printf("%d\n", INTEGER(nrm)[0]);
// 	INTEGER_POINTER(getListElement(bidirOptimParams, "nrm"))[LENGTH(getListElement(bidirOptimParams, "nrm"))] = INTEGER(nrm)[0];
// 	for(i=0; i<LENGTH(getListElement(bidirOptimParams, "nrm")); i++) {
// 		printf("%d\n", INTEGER(getListElement(bidirOptimParams, "nrm"))[i]);
// 	}

//UNPROTECT(2) ;
//}
    UNPROTECT(2) ;

    for(i=0; i<this->K; i++)
    {
        for(j=0; j<this->K; j++)
        {
            this->updateNumerator[i][j] = 0;
            this->updateDenominator[i][j] = 0;
        }
    }

//this->callRsolnp(bidirOptimParams);
//this->finalize();
    if(DEBUG)
    {
        Rprintf("Optimization finished.\n");
    }
}


SEXP TransitionMatrix::callRsolnp(SEXP pars)
{

    SEXP oldpars = getListElement(pars, "x0");

    int i;
    for(i=0; i<LENGTH(oldpars); i++)
    {
        Rprintf("%f ", REAL(oldpars)[i]);
    }
    Rprintf("\n");

    SEXP call = PROTECT( lang2( install( "c2solnp"), pars ) ) ;
    SEXP res = PROTECT( eval( call, R_GlobalEnv ) ) ;

    SEXP newpars = getListElement(res, "pars");

    int convergence = (int)REAL(getListElement(res, "convergence"))[0];
    if(convergence == 0)
    {
        SET_ELEMENT(pars, 0, newpars);
        for(i=0; i<LENGTH(VECTOR_ELT(pars, 0)); i++)
        {
            Rprintf("%f ", REAL(VECTOR_ELT(pars, 0))[i]);
        }
        Rprintf("\n");
    }
    else
    {
        SEXP warn = PROTECT( lang2( install( "warning"), mkString("Rsolnp did not converge. Using old parameter estimates!\n") ) ) ;
        SEXP exwarn = PROTECT( eval( warn, R_GlobalEnv ) ) ;
        UNPROTECT(2);
    }

    UNPROTECT(2) ;

    return res ;
}


// slight numerical instabilities
void TransitionMatrix::finalize()
{
    int i,j;
    for (i=0; i<this->K; i++)
    {
        double rowSum=0;
        for (j=0; j<this->K; j++)
        {
            rowSum += this->A[i][j];
        }
        for (j=0; j<this->K; j++)
        {
            this->A[i][j] /= rowSum;
        }
    }
}


int TransitionMatrix::getK()
{
    return this->K;
}
