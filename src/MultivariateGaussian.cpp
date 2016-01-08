#include "MultivariateGaussian.h"
#include <cmath>
#define OBSERVATIONPOS(index, start) (start[index])

MultivariateGaussian::MultivariateGaussian(
ParamContainerEmissions *emissionParams)
{

    this->emissionParams = emissionParams;

    int mem = sizeof(double) * this->emissionParams->getD() * 4
        + sizeof(double*) * this->emissionParams->getD() * 2;

    if (DEBUG_MEMORY)
    {
        printf("new->MultivariateGaussian; (%d bytes) \n", mem);
    }

    this->updateNumeratorMU = (double*) malloc(
        sizeof(double) * this->emissionParams->getD());
    this->updateDenominatorMU = (double*) malloc(
        sizeof(double) * this->emissionParams->getD());
    this->updateNumeratorSIGMA = (double**) malloc(
        sizeof(double*) * this->emissionParams->getD());
    this->updateDenominatorSIGMA = (double**) malloc(
        sizeof(double*) * this->emissionParams->getD());

    int d1, d2;
    for (d1 = 0; d1 < this->emissionParams->getD(); d1++)
    {
        this->updateNumeratorMU[d1] = 0;
        this->updateDenominatorMU[d1] = 0;
        this->updateNumeratorSIGMA[d1] = (double*) malloc(
            sizeof(double) * this->emissionParams->getD());
        this->updateDenominatorSIGMA[d1] = (double*) malloc(
            sizeof(double) * this->emissionParams->getD());
        for (d2 = 0; d2 < this->emissionParams->getD(); d2++)
        {
            this->updateNumeratorSIGMA[d1][d2] = 0;
            this->updateDenominatorSIGMA[d1][d2] = 0;
        }
    }

}


MultivariateGaussian::~MultivariateGaussian()
{

    int d1;
    free(this->updateNumeratorMU);
    free(this->updateDenominatorMU);

    for (d1 = 0; d1 < this->emissionParams->getD(); d1++)
    {
        free(this->updateNumeratorSIGMA[d1]);
        free(this->updateDenominatorSIGMA[d1]);
    }
    free(this->updateNumeratorSIGMA);
    free(this->updateDenominatorSIGMA);
    int mem = sizeof(double) * this->emissionParams->getD() * 4
        + sizeof(double*) * this->emissionParams->getD() * 2;
    if (DEBUG_MEMORY)
    {
        printf("delete->MultivariateGaussian; (%d bytes) ", mem);
    }
    delete this->emissionParams;
}


double MultivariateGaussian::calcEmissionProbability(double *obs, int isna, int currN)
{
//Rprintf("MG\n");
//std::cout << "Using MultivariateGaussian\n";
    int i, j, D;
    double exponent = 0;
    double first_term = 0;
    int obs_i = 0;
    int obs_j = 0;
    first_term = pow(2.5066282746310002, this->emissionParams->getD())
        * sqrt(this->emissionParams->getGaussianDET());
    double density = 1;
    int* myStart = this->emissionParams->getStart();
    D = this->emissionParams->getD();
    double** myMU = this->emissionParams->getGaussianMU();
    double** myInvSigma = this->emissionParams->getGaussianINVSIGMA();
    if (isna == 0)
    {
        for (i = 0; i < D; i++)
        {
            obs_i = myStart[i];
//printf("obs_i=%d , obs=%f dim=%d\n", obs_i, obs[obs_i], this->emissionParams->getD());
            if (obs[obs_i] != obs[obs_i])
            {
                isna = 1;
                break;
            }
            for (j = 0; j < D; j++)
            {

                obs_j = myStart[j];
//	printf("j=%d, obs_j=%d , obs=%f dim=%d\n", j, obs_j, obs[obs_j], this->emissionParams->getD());
                if (obs[obs_j] != obs[obs_j])
                {
                    isna = 1;
                    break;
                }
                exponent += (obs[obs_i]
                    - myMU[i][0])
                    * myInvSigma[i][j]
                    * (obs[obs_j]
                    - myMU[j][0]);
            }
        }
        density = exp(-0.5 * exponent) / first_term;
    }
    if (density > 1e20)
    {
        error("Ill-conditioned covariance matrix!\n");
    }
//	Rprintf("dmg=%f\n", density);
//	density =  density*pow(2.5066282746310002,  this->emissionParams->getLogCovPrior() );

    if (density < 1e-100)
    {
        density = 1e-100;
    }

    return density;
}


void MultivariateGaussian::updateAuxiliaries(double*** observations,
double** gamma, double* Pk, int* T, int n, int i, int** isNaN)
{
    int t, d, l;
    double numer, denom;
    int obs_d = 0;
//Rprintf("Gauss-aux\n");
    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        obs_d = OBSERVATIONPOS(d, this->emissionParams->getStart());
        numer = 0.0;
        denom = 0.0;
        for (t = 0; t < T[n]; t++)
        {
            if (isNaN[n][t] == 0)
            {
                numer = numer + gamma[t][i] * observations[n][t][obs_d];
                denom = denom + gamma[t][i];
            }
        }
        this->updateNumeratorMU[d] += 1 / Pk[n] * numer;
        this->updateDenominatorMU[d] += 1 / Pk[n] * denom;
//Rprintf("%f/%f\n", this->updateNumeratorMU[d], this->updateDenominatorMU[d]);
    }

    for (t = 0; t < T[n]; t++)
    {
        this->emissionParams->setGammaAux(gamma[t][i], n, t);
    }
}


void MultivariateGaussian::updateAuxiliariesCoupled(double*** observations,
double** gamma, double* Pk, int* T, int n, int i, int statecouple,
int** isNaN)
{
    int t, d, l, obs_d;
    double numer, denom;

    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        numer = 0.0;
        denom = 0.0;
        obs_d = OBSERVATIONPOS(d, this->emissionParams->getStart());
        for (t = 0; t < T[n]; t++)
        {
            if (isNaN[n][t] == 0)
            {
                numer = numer
                    + (gamma[t][i] + gamma[t][statecouple])
                    * observations[n][t][obs_d];
                denom = denom + (gamma[t][i] + gamma[t][statecouple]);
            }
        }
        this->updateNumeratorMU[d] += 1 / Pk[n] * numer;
        this->updateDenominatorMU[d] += 1 / Pk[n] * denom;
    }

    for (t = 0; t < T[n]; t++)
    {
        this->emissionParams->setGammaAux((gamma[t][i] + gamma[t][statecouple]),
            n, t);
    }
}


void MultivariateGaussian::updateAuxiliariesCoupledRevop(double*** observations,
double** gamma, double* Pk, int* T, int n, int i, int statecouple,
int* state2flag, int* revop, int** isNaN)
{
    int t, d, l, obs_d;
    double numer, denom;

//printf(" GAUSSIAN s=%d, c=%d\n", i, statecouple);
    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        numer = 0.0;
        denom = 0.0;
        obs_d = OBSERVATIONPOS(d, this->emissionParams->getStart());
//printf("GAUSSIAN d= %d obs_d=%d \n",d, obs_d);
        for (t = 0; t < T[n]; t++)
        {
            if (isNaN[n][t] == 0)
            {

                if (state2flag[statecouple] == 1)
                {
                    numer = numer + gamma[t][i] * observations[n][t][obs_d]
                        + gamma[t][statecouple]
                                                  //// d
                        * observations[n][t][revop[obs_d]];
                }
                else
                {
                    numer = numer
                        + gamma[t][i] * observations[n][t][revop[obs_d]]
                        + gamma[t][statecouple] * observations[n][t][obs_d];
                }
                denom = denom + (gamma[t][i] + gamma[t][statecouple]);
            }
        }
        this->updateNumeratorMU[d] += 1 / Pk[n] * numer;
        this->updateDenominatorMU[d] += 1 / Pk[n] * denom;

    }

    for (t = 0; t < T[n]; t++)
    {
        this->emissionParams->setGammaAux(gamma[t][i], n, t);
    }
}


void MultivariateGaussian::updateCoupledRevop(double ***observations,
double* Pk, int statecouple, int* state2flag, int* revop,
double** revGammaAux, int** isNaN, SEXP emissionPrior, int currN, int ncores)
{

    int skipcounter = 0;
    int d;
    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        this->emissionParams->setGaussianMUelement(
            this->updateNumeratorMU[d] / this->updateDenominatorMU[d], d);
        this->updateNumeratorMU[d] = 0;
        this->updateDenominatorMU[d] = 0;
//printf("%f ", this->emissionParams->getGaussianMU()[d][0]);
    }
/*for(d=0; d<6;d++) {
Rprintf("%d ", revop[d]);
}
Rprintf("\n");*/
    double numer, denom;
    double** myTempNumer = (double**) malloc(sizeof(double*) * this->emissionParams->getD());
    double** myTempDenom = (double**) malloc(sizeof(double*) * this->emissionParams->getD());

    int i, j;
    for(i = 0; i < this->emissionParams->getD(); i++)
    {
        myTempNumer[i] = (double*) malloc(sizeof(double) * this->emissionParams->getD());
        myTempDenom[i] = (double*) malloc(sizeof(double) * this->emissionParams->getD());
        for (j = 0; j < this->emissionParams->getD(); j++)
        {
            myTempNumer[i][j] = 0;
            myTempDenom[i][j] = 0;
        }
    }

    int n;

    int lower_n = 0;
    int upper_n = this->emissionParams->getNsample();
    if(currN != -1)
    {
        lower_n = currN;
        upper_n = currN+1;
    }

    if(this->emissionParams->getD() < ncores)
    {
        ncores = this->emissionParams->getD();
    }
    int *myDimBuckets = (int*)malloc(sizeof(int)*ncores+1);
    for(i=0; i<=ncores; i++)
    {
        myDimBuckets[i] = 0;
    }
    int currbucket = 1;
    for(i=0; i<this->emissionParams->getD(); i++)
    {
        myDimBuckets[currbucket] = myDimBuckets[currbucket]+1;
        if(currbucket == ncores)
        {
            currbucket = 0;
        }
        currbucket++;
    }
    for(i=1; i<ncores+1; i++)
    {
        myDimBuckets[i] = myDimBuckets[i]+myDimBuckets[i-1];
    }

    if(this->emissionParams->getUpdateCov())
    {
//Rprintf("Updating covariances\n");
        int* myStart = this->emissionParams->getStart();
        int myD = this->emissionParams->getD();
        int* myT = this->emissionParams->getT();
        double** myGammaAux = this->emissionParams->getGammaAux();
        double** myMU = this->emissionParams->getGaussianMU();
        for(n = lower_n; n < upper_n; n++)
        {
            int k;
#pragma omp parallel for
            for(k=1; k<ncores+1; k++)
            {
                int d1;
                for(d1=myDimBuckets[k-1]; d1<myDimBuckets[k]; d1++)
                {
                    int obs_d1 = myStart[d1];
                    int d2;
                    for (d2 = d1; d2 < myD; d2++)
                    {
                        int obs_d2 = myStart[d2];
                        int t, l;
                        for (t = 0; t < myT[n]; t++)
                        {
                            if (isNaN[n][t] == 0)
                            {
                                if (state2flag[statecouple] == 1)
                                {
                                    myTempNumer[d1][d2] += myGammaAux[n][t] * (observations[n][t][obs_d1] - myMU[d1][0]) * (observations[n][t][obs_d2] - myMU[d2][0])
                                        + revGammaAux[n][t] * (observations[n][t][revop[obs_d1]] - myMU[d1][0]) * (observations[n][t][revop[obs_d2]] - myMU[d2][0]);
                                }
                                else
                                {
                                    myTempNumer[d1][d2] += myGammaAux[n][t] * (observations[n][t][revop[obs_d1]] - myMU[d1][0]) * (observations[n][t][revop[obs_d2]] - myMU[d2][0])
                                        + revGammaAux[n][t] * (observations[n][t][obs_d1] - myMU[d1][0]) * (observations[n][t][obs_d2] - myMU[d2][0]);
                                }
                                myTempDenom[d1][d2] += myGammaAux[n][t] + revGammaAux[n][t];
                            }

                        }
                        this->updateNumeratorSIGMA[d1][d2] += 1.0 / Pk[n] * myTempNumer[d1][d2];
                        this->updateDenominatorSIGMA[d1][d2] += 1.0 / Pk[n] * myTempDenom[d1][d2];
                        if (d1 != d2)
                        {
                            this->updateNumeratorSIGMA[d2][d1] = this->updateNumeratorSIGMA[d1][d2];
                            this->updateDenominatorSIGMA[d2][d1] = this->updateDenominatorSIGMA[d1][d2];
                        }
//	Rprintf("myTempNumer[%d][%d]=%f\n", d1,d2, myTempNumer[d1][d2]);

                    }
                }
            }
        }

        for (i = 0; i < this->emissionParams->getD(); i++)
        {
            for (j = 0; j < this->emissionParams->getD(); j++)
            {
                this->emissionParams->setGaussianSIGMAelement(this->updateNumeratorSIGMA[i][j]/this->updateDenominatorSIGMA[i][j], i, j);
                this->emissionParams->setGaussianINVSIGMAelement(this->updateNumeratorSIGMA[i][j]/this->updateDenominatorSIGMA[i][j], i, j);
//	Rprintf("%f (%f) (%f)", this->emissionParams->getGaussianSIGMA()[i][j], this->updateNumeratorSIGMA[i][j], this->updateDenominatorSIGMA[i][j]);
                if(this->emissionParams->getSharedCov() == 0)
                {
                    this->updateNumeratorSIGMA[i][j] = 0;
                    this->updateDenominatorSIGMA[i][j] = 0;
                }
            }
//	Rprintf("\n");
        }
//Rprintf("\n\n");
        inverse(this->emissionParams->getGaussianINVSIGMA(), this->emissionParams->getD());
        this->emissionParams->setGaussianDET(matrixDet(this->emissionParams->getGaussianSIGMA(), this->emissionParams->getD()));

    }

    for (i = 0; i < this->emissionParams->getD(); i++)
    {
        free(myTempNumer[i]);
        free(myTempDenom[i]);
    }
    free(myTempNumer);
    free(myTempDenom);
    free(myDimBuckets);

}


double MultivariateGaussian::Prior(SEXP hyperparams)
{
    int i, j;
    for (i = 0; i < this->emissionParams->getD(); i++)
    {
        for (j = 0; j < this->emissionParams->getD(); j++)
        {
            REAL(getListElement(hyperparams, "cov"))[i
                + j * this->emissionParams->getD()] =
                this->emissionParams->getGaussianSIGMA()[i][j];
// Rprintf("%f ", this->emissionParams->getGaussianSIGMA()[i][j]);
        }
//Rprintf("\n");
    }

// call solnp from R for optimization
//SEXP call = PROTECT(lang2(install("calldiwish"), hyperparams));
    SEXP call = PROTECT(lang2(getListElement(hyperparams, "calldiwish"), hyperparams));

    SEXP res = PROTECT(eval(call, R_GlobalEnv));
    double out = REAL(res)[0];
//Rprintf("%f\n", out);
    UNPROTECT(2);
    return out;

}


void MultivariateGaussian::update(double ***observations, double* Pk,
int** isNaN, SEXP emissionPrior, int currN, int ncores)
{

    int d;
    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        this->emissionParams->setGaussianMUelement(
            this->updateNumeratorMU[d] / this->updateDenominatorMU[d], d);
        this->updateNumeratorMU[d] = 0;
        this->updateDenominatorMU[d] = 0;
//printf("%f ", this->emissionParams->getGaussianMU()[d][0]);
    }

    double numer, denom;
    double** myTempNumer = (double**) malloc(sizeof(double*) * this->emissionParams->getD());
    double** myTempDenom = (double**) malloc(sizeof(double*) * this->emissionParams->getD());

    int i, j;
    for(i = 0; i < this->emissionParams->getD(); i++)
    {
        myTempNumer[i] = (double*) malloc(sizeof(double) * this->emissionParams->getD());
        myTempDenom[i] = (double*) malloc(sizeof(double) * this->emissionParams->getD());
        for (j = 0; j < this->emissionParams->getD(); j++)
        {
            myTempNumer[i][j] = 0;
            myTempDenom[i][j] = 0;
        }
    }

    int n;

    int lower_n = 0;
    int upper_n = this->emissionParams->getNsample();
    if(currN != -1)
    {
        lower_n = currN;
        upper_n = currN+1;
    }

    if(this->emissionParams->getD() < ncores)
    {
        ncores = this->emissionParams->getD();
    }
    int *myDimBuckets = (int*)malloc(sizeof(int)*ncores+1);
    for(i=0; i<=ncores; i++)
    {
        myDimBuckets[i] = 0;
    }
    int currbucket = 1;
    for(i=0; i<this->emissionParams->getD(); i++)
    {
        myDimBuckets[currbucket] = myDimBuckets[currbucket]+1;
        if(currbucket == ncores)
        {
            currbucket = 0;
        }
        currbucket++;
    }
    for(i=1; i<ncores+1; i++)
    {
        myDimBuckets[i] = myDimBuckets[i]+myDimBuckets[i-1];
    }

    if(this->emissionParams->getUpdateCov())
    {
//Rprintf("Updating covariances\n");
        int* myStart = this->emissionParams->getStart();
        int myD = this->emissionParams->getD();
        int* myT = this->emissionParams->getT();
        double** myGammaAux = this->emissionParams->getGammaAux();
        double** myMU = this->emissionParams->getGaussianMU();

        for(n = lower_n; n < upper_n; n++)
        {
            int k;
#pragma omp parallel for
            for(k=1; k<ncores+1; k++)
            {
                int d1;
                for(d1=myDimBuckets[k-1]; d1<myDimBuckets[k]; d1++)
                {
                    int obs_d1 = myStart[d1];
                    int d2;
                    for (d2 = d1; d2 < myD; d2++)
                    {
                        int obs_d2 = myStart[d2];
                        int t, l;
                        for (t = 0; t < myT[n]; t++)
                        {
                            if (isNaN[n][t] == 0)
                            {
                                myTempNumer[d1][d2] += myGammaAux[n][t] * (observations[n][t][obs_d1] - myMU[d1][0]) * (observations[n][t][obs_d2] - myMU[d2][0]);
                                myTempDenom[d1][d2] += myGammaAux[n][t];
                            }

                        }

                        this->updateNumeratorSIGMA[d1][d2] += 1.0 / Pk[n] * myTempNumer[d1][d2];
                        this->updateDenominatorSIGMA[d1][d2] += 1.0 / Pk[n] * myTempDenom[d1][d2];
                        if (d1 != d2)
                        {
                            this->updateNumeratorSIGMA[d2][d1] = this->updateNumeratorSIGMA[d1][d2];
                            this->updateDenominatorSIGMA[d2][d1] = this->updateDenominatorSIGMA[d1][d2];
                        }
                    }
                }
            }
        }

        for (i = 0; i < this->emissionParams->getD(); i++)
        {
            for (j = 0; j < this->emissionParams->getD(); j++)
            {
                this->emissionParams->setGaussianSIGMAelement(this->updateNumeratorSIGMA[i][j]/this->updateDenominatorSIGMA[i][j], i, j);
                this->emissionParams->setGaussianINVSIGMAelement(this->updateNumeratorSIGMA[i][j]/this->updateDenominatorSIGMA[i][j], i, j);
                if(this->emissionParams->getSharedCov() == 0)
                {
                    this->updateNumeratorSIGMA[i][j] = 0;
                    this->updateDenominatorSIGMA[i][j] = 0;
                }
            }
        }

        inverse(this->emissionParams->getGaussianINVSIGMA(), this->emissionParams->getD());
        this->emissionParams->setGaussianDET(matrixDet(this->emissionParams->getGaussianSIGMA(), this->emissionParams->getD()));

    }

    for (i = 0; i < this->emissionParams->getD(); i++)
    {
        free(myTempNumer[i]);
        free(myTempDenom[i]);
    }
    free(myTempNumer);
    free(myTempDenom);
    free(myDimBuckets);

}


double** MultivariateGaussian::getUpdateNumeratorSigma()
{
    return this->updateNumeratorSIGMA;
}


double** MultivariateGaussian::getUpdateDenominatorSigma()
{
    return this->updateDenominatorSIGMA;
}


void MultivariateGaussian::computeShared(EmissionFunction** myEmissions, int nStates)
{
//Rprintf("sC=%d\n", this->emissionParams->getSharedCov());
    if(this->emissionParams->getSharedCov() == 1)
    {

        double** myTempNumer = (double**) malloc(
            sizeof(double*) * this->emissionParams->getD());
        double** myTempDenom = (double**) malloc(
            sizeof(double*) * this->emissionParams->getD());

        int d1, d2;
        for (d1 = 0; d1 < this->emissionParams->getD(); d1++)
        {
            myTempNumer[d1] = (double*) malloc(
                sizeof(double) * this->emissionParams->getD());
            myTempDenom[d1] = (double*) malloc(
                sizeof(double) * this->emissionParams->getD());
            for (d2 = 0; d2 < this->emissionParams->getD(); d2++)
            {
                myTempNumer[d1][d2] = 0;
                myTempDenom[d1][d2] = 0;
            }
        }
//Rprintf("1\n");
        int k;
        for(k=0; k<nStates; k++)
        {
//		printf("k=%d\n",k);
            for (d1 = 0; d1 < this->emissionParams->getD(); d1++)
            {
                for (d2 = 0; d2 < this->emissionParams->getD(); d2++)
                {
//						printf("whichone=%d, w2=%d\n", myEmissions[k]->getParameter()->getWhichOne(), this->getParameter()->getWhichOne());
                    if(myEmissions[k]->getParameter()->getWhichOne() == 1)
                    {
//							printf("k=%d, d1=%d, d2=%d, myD=%d\n", k, d1, d2, this->emissionParams->getD());
                        myTempNumer[d1][d2] += myEmissions[k]->updateNumeratorSIGMA[d1][d2];
                        myTempDenom[d1][d2] += myEmissions[k]->updateDenominatorSIGMA[d1][d2];
                    }

                }
            }
        }

        for (d1 = 0; d1 < this->emissionParams->getD(); d1++)
        {
            for (d2 = 0; d2 < this->emissionParams->getD(); d2++)
            {
                this->emissionParams->setGaussianSIGMAelement((myTempNumer[d1][d2] / myTempDenom[d1][d2]), d1, d2);
                this->emissionParams->setGaussianINVSIGMAelement((myTempNumer[d1][d2] / myTempDenom[d1][d2]), d1, d2);
//					Rprintf("%f ", this->emissionParams->getGaussianSIGMA()[d1][d2]);
            }
//				Rprintf("\n");
        }
//		Rprintf("\n");
        inverse(this->emissionParams->getGaussianINVSIGMA(), this->emissionParams->getD());
        this->emissionParams->setGaussianDET(matrixDet(this->emissionParams->getGaussianSIGMA(), this->emissionParams->getD()));
        for (d1 = 0; d1 < this->emissionParams->getD(); d1++)
        {
            free(myTempNumer[d1]);
            free(myTempDenom[d1]);
        }
        free(myTempNumer);
        free(myTempDenom);
    }
//Rprintf("2\n");
}


void MultivariateGaussian::resetShared()
{
//	Rprintf("fuck\n");
    int d1,d2;
    for (d1 = 0; d1 < this->emissionParams->getD(); d1++)
    {
        for (d2 = 0; d2 < this->emissionParams->getD(); d2++)
        {
            this->updateNumeratorSIGMA[d1][d2] = 0;
            this->updateDenominatorSIGMA[d1][d2] = 0;
        }
    }
//Rprintf("fuck2\n");

}


void MultivariateGaussian::setParsToTwin(EmissionFunction* myTwinEmission, int currN, double*** observations)
{

    if(this->emissionParams->getUpdateCov())
    {
        int d;
        for (d = 0; d < this->emissionParams->getD(); d++)
        {
            this->emissionParams->setGaussianMUelement(myTwinEmission->getParameter()->getGaussianMU()[d][0], d);
            this->updateNumeratorMU[d] = 0;
            this->updateDenominatorMU[d] = 0;
        }

        int d1, d2;
        for (d1 = 0; d1 < this->emissionParams->getD(); d1++)
        {
            for (d2 = 0; d2 < this->emissionParams->getD(); d2++)
            {
                this->emissionParams->setGaussianSIGMAelement(myTwinEmission->getParameter()->getGaussianSIGMA()[d1][d2], d1, d2);
                this->emissionParams->setGaussianINVSIGMAelement(myTwinEmission->getParameter()->getGaussianINVSIGMA()[d1][d2], d1, d2);
                if(this->emissionParams->getSharedCov() == 0)
                {
                    this->updateNumeratorSIGMA[d1][d2] = 0;
                    this->updateDenominatorSIGMA[d1][d2] = 0;
                }
            }
        }
        this->emissionParams->setGaussianDET(myTwinEmission->getParameter()->getGaussianDET());

    }

}
