
#include "Multinomial.h"
#include <cmath>
#include <list>
using namespace std;

#define OBSERVATIONPOS(index, start) (start[index])

Multinomial::Multinomial(ParamContainerEmissions *emissionParams)
{

    this->emissionParams = emissionParams;

    int mem = sizeof(double) * this->emissionParams->getD() * 4
        + sizeof(double*) * this->emissionParams->getD() * 2;
    if (DEBUG_MEMORY)
    {
        printf("new->Multinomial; (%d bytes) \n", mem);
    }

    this->updateNumeratorMP = (double*) malloc(
        sizeof(double) * this->emissionParams->getD());
    this->updateDenominatorMP = (double*) malloc(
        sizeof(double) * this->emissionParams->getD());

    int d1, d2;
    for (d1 = 0; d1 < this->emissionParams->getD(); d1++)
    {
        this->updateNumeratorMP[d1] = 0;
        this->updateDenominatorMP[d1] = 0;
    }
}


Multinomial::~Multinomial()
{

    int d1;
    free(this->updateNumeratorMP);
    free(this->updateDenominatorMP);

    int mem = sizeof(double) * this->emissionParams->getD() * 4
        + sizeof(double*) * this->emissionParams->getD() * 2;
    if (DEBUG_MEMORY)
    {
        printf("delete->Multinomial; (%d bytes) ", mem);
    }
    delete this->emissionParams;
//delete this->emissionParams;
}


double Multinomial::calcEmissionProbability(double *obs, int isna, int currN)
{
//printf("Multinomial\n");
    int myN = 0;                                  //this->emissionParams->getN();
    int myD = this->emissionParams->getD();
    double mult_coeff = 0;
    double numer = 0;
    double denom = 0;
    double probability = 1;

    int* revComp = this->emissionParams->getReverseComplementary();
    int myDir = this->emissionParams->getStateDir();

//Rprintf("%d\n", myDir);
    int i, d, obs_d;
    for(d=0; d<myD; d++)
    {
        if(myDir == -1)
        {
            obs_d = OBSERVATIONPOS(revComp[d], this->emissionParams->getStart());
        }
        else
        {
            obs_d = OBSERVATIONPOS(d, this->emissionParams->getStart());
        }

//		Rprintf("D=%d, obs_d=%d, d=%d, obs=%f\n", myD, obs_d, d, obs[obs_d]);
        myN += obs[obs_d];
    }
//	printf("myN=%d\n", myN);
    if (isna == 0 && myN >= 1)
    {

        for(i=1; i<=myN; i++)
        {
            numer += log(i);
        }

        for(d=0; d<myD; d++)
        {
            if(myDir == -1)
            {
                obs_d = OBSERVATIONPOS(revComp[d], this->emissionParams->getStart());
            }
            else
            {
                obs_d = OBSERVATIONPOS(d, this->emissionParams->getStart());
            }

            for(i=1; i<=obs[obs_d]; i++)
            {
                denom += log(i);
            }
        }

        double* myP = this->emissionParams->getMultinomialP();
        double logProbSum = 0;
        for(d=0; d<myD; d++)
        {
            if(myDir == -1)
            {
                obs_d = OBSERVATIONPOS(revComp[d], this->emissionParams->getStart());
            }
            else
            {
                obs_d = OBSERVATIONPOS(d, this->emissionParams->getStart());
            }

            logProbSum += obs[obs_d]*log(myP[d]+1e-300);
//	printf("obs=%f, p=%f, log(p)=%f\n", obs[obs_d], myP[d], log(myP[d]+1e-300));
//	for(i=0; i<obs[obs_d]; i++) {
//		 logProbSum += log(myP[(int)obs[obs_d]]);
//	}
        }
//	Rprintf("%f, %f, %f, %f\n", obs[0], obs[1], exp((numer - denom) + logProbSum), (numer - denom) + logProbSum);

        probability = exp((numer - denom) + logProbSum);
        if ((probability != probability) | (probability > 1e20))
        {
            warning("Multinomial emission probability calculation is instable.");
//Rprintf("%f, %f, %f, %f\n", numer, denom, logProbSum, (numer - denom) + logProbSum);
        }
    }
    else
    {
//Rprintf("NA here!\n");
    }

//	Rprintf("dmn=%f\n", probability);
    if (probability < 1e-100)
    {
//Rprintf("no no no\n");
        probability = 1e-100;
    }
    return probability;
}


void Multinomial::updateAuxiliaries(double*** observations, double** gamma,
double* Pk, int* T, int n, int i, int** isNaN)
{
    int t, d, l, obs_d, obs_d_revcomp;
    double numer, denom;
//	Rprintf("Mult-aux\n");
    int* revComp = this->emissionParams->getReverseComplementary();

    double* mySums = (double*)malloc(sizeof(double)*T[n]);
    for (t = 0; t < T[n]; t++)
    {
        mySums[t] = 0;
        for (d = 0; d < this->emissionParams->getD(); d++)
        {
            obs_d = OBSERVATIONPOS(d, this->emissionParams->getStart());
            obs_d_revcomp = OBSERVATIONPOS(revComp[d], this->emissionParams->getStart());

            mySums[t] += observations[n][t][obs_d]+observations[n][t][obs_d_revcomp];
        }
    }

    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        obs_d = OBSERVATIONPOS(d, this->emissionParams->getStart());
        obs_d_revcomp = OBSERVATIONPOS(revComp[d], this->emissionParams->getStart());

        numer = 0.0;
        denom = 0.0;
        for (t = 0; t < T[n]; t++)
        {

            if (isNaN[n][t] == 0)
            {
//numer = numer + gamma[t][i] * observations[n][t][d]; //indikator
                                                  //IndicatorFct(observations[n][t][obs_d]);
                numer = numer + gamma[t][i] * (observations[n][t][obs_d]+observations[n][t][obs_d_revcomp]);
                denom = denom + gamma[t][i]*mySums[t];
            }
//printf("obs[%d][%d]=%f, gamma=%f\n", t, obs_d, observations[n][t][obs_d],gamma[t][i]);
        }
        this->updateNumeratorMP[d] += 1 / Pk[n] * numer;
        this->updateDenominatorMP[d] += 1 / Pk[n] * denom;

    }
    free(mySums);
//printf("\n Multinomial type %f/%f\n", this->updateDenominatorMP[0], this->updateNumeratorMP[0] );
}


void Multinomial::updateAuxiliariesCoupled(double*** observations, double** gamma,
double* Pk, int* T, int n, int i, int statecouple, int** isNaN)
{
    error("Should bot be here!");

}


void Multinomial::updateAuxiliariesCoupledRevop(double*** observations,
double** gamma, double* Pk, int* T, int n, int i, int statecouple,
int* state2flag, int* revop, int** isNaN)
{
    int t, d, l, obs_d,obs_d_revcomp;
    double numer, denom;
    int* revComp = this->emissionParams->getReverseComplementary();

    double* mySums = (double*)malloc(sizeof(double)*T[n]);
    for (t = 0; t < T[n]; t++)
    {
        mySums[t] = 0;
        for (d = 0; d < this->emissionParams->getD(); d++)
        {
            obs_d = OBSERVATIONPOS(d, this->emissionParams->getStart());
            obs_d_revcomp = OBSERVATIONPOS(revComp[d], this->emissionParams->getStart());
//			printf("obs_d=%d, revobs_d=%d, revd=%d\n", obs_d, obs_d_revcomp, revComp[d]);

            if (state2flag[statecouple] == 1)
            {
                mySums[t] += (gamma[t][i]*observations[n][t][obs_d]) + (gamma[t][statecouple]*observations[n][t][revop[obs_d_revcomp]]);

            }
            else
            {
                mySums[t] += (gamma[t][i]*observations[n][t][revop[obs_d_revcomp]]) + (gamma[t][statecouple]*observations[n][t][obs_d]);
            }

        }
    }

    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        obs_d = OBSERVATIONPOS(d, this->emissionParams->getStart());
        obs_d_revcomp = OBSERVATIONPOS(revComp[d], this->emissionParams->getStart());
        numer = 0.0;
        denom = 0.0;
        for (t = 0; t < T[n]; t++)
        {
            if (isNaN[n][t] == 0)
            {
                if (state2flag[statecouple] == 1)
                {
                    numer = numer + (gamma[t][i]*observations[n][t][obs_d])
                        + (gamma[t][statecouple]*observations[n][t][revop[obs_d_revcomp]]);

                }
                else
                {
                    numer = numer
                        + (gamma[t][i]*observations[n][t][revop[obs_d_revcomp]])
                        + (gamma[t][statecouple]*observations[n][t][obs_d]);
                }
                denom = denom + mySums[t];
            }
        }
        this->updateNumeratorMP[d] += 1 / Pk[n] * numer;
        this->updateDenominatorMP[d] += 1 / Pk[n] * denom;
    }
    free(mySums);
}


void Multinomial::updateCoupledRevop(double ***observations, double* Pk,
int statecouple, int* state2flag, int* revop, double** revGammaAux,
int** isNaN, SEXP emissionPrior, int currN, int ncores)
{

    int d;
    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        this->emissionParams->setMultinomialPelement(
            this->updateNumeratorMP[d] / this->updateDenominatorMP[d], d);

        this->updateNumeratorMP[d] = 0;
        this->updateDenominatorMP[d] = 0;

    }

}


// was mach ich aus dem prior
double Multinomial::Prior(SEXP hyperparams)
{
    return 0;
}


void Multinomial::update(double ***observations, double* Pk, int** isNaN,
SEXP emissionPrior, int currN, int ncores)
{
    int d;
    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        this->emissionParams->setMultinomialPelement(
            this->updateNumeratorMP[d] / this->updateDenominatorMP[d], d);

        this->updateNumeratorMP[d] = 0;
        this->updateDenominatorMP[d] = 0;

    }
}


void Multinomial::setParsToTwin(EmissionFunction* myTwinEmission, int currN, double*** observations)
{
    int d;
    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        this->emissionParams->setMultinomialPelement(myTwinEmission->getParameter()->getMultinomialP()[d], d);

        this->updateNumeratorMP[d] = 0;
        this->updateDenominatorMP[d] = 0;

    }
}
