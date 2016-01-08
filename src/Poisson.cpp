
#include "Poisson.h"
#include <cmath>
#include <list>
using namespace std;

#define OBSERVATIONPOS(index, start) (start[index])

Poisson::Poisson(ParamContainerEmissions *emissionParams)
{

    this->emissionParams = emissionParams;

    int mem = sizeof(double) * this->emissionParams->getD() * 4
        + sizeof(double*) * this->emissionParams->getD() * 2;
    if (DEBUG_MEMORY)
    {
        printf("new->Poisson; (%d bytes) \n", mem);
    }

    this->updateNumeratorLambda = (double*) malloc(
        sizeof(double) * this->emissionParams->getD());
    this->updateDenominatorLambda = (double*) malloc(
        sizeof(double) * this->emissionParams->getD());

    int d1, d2;
    for (d1 = 0; d1 < this->emissionParams->getD(); d1++)
    {
        this->updateNumeratorLambda[d1] = 0;
        this->updateDenominatorLambda[d1] = 0;
    }
}


Poisson::~Poisson()
{

    int d1;
    free(this->updateNumeratorLambda);
    free(this->updateDenominatorLambda);

    int mem = sizeof(double) * this->emissionParams->getD() * 4
        + sizeof(double*) * this->emissionParams->getD() * 2;
    if (DEBUG_MEMORY)
    {
        printf("delete->Poisson; (%d bytes) ", mem);
    }
//delete this->emissionParams;
}


double Poisson::getP()
{
    return this->emissionParams->getPoissonLambda();
}


double Poisson::calcEmissionProbability(double *obs, int isna, int currN)
{

    int i, j, obs_i;
    double probability = 1;
    double myLambda = this->emissionParams->getPoissonLambda();
    if (isna == 0)
    {
        for(i = 0; i < this->emissionParams->getD(); i++)
        {
            obs_i = OBSERVATIONPOS(i, this->emissionParams->getStart());

            if (obs[obs_i] != obs[obs_i])         //check NAN
            {
                isna = 1;
                break;
            }
            else
            {
                double currK = obs[obs_i];
                for(j=1; j<=currK; j++)
                {
                    probability = probability*(myLambda/j);
                }
                probability = probability
                    * exp(-myLambda);

            }
        }

    }

    if (probability < 1e-100)
    {
        probability = 1e-100;
    }
    return probability;
}


void Poisson::updateAuxiliaries(double*** observations, double** gamma,
double* Pk, int* T, int n, int i, int** isNaN)
{
    int t, d, l, obs_d;
    double numer, denom;
    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        obs_d = OBSERVATIONPOS(d, this->emissionParams->getStart());

        numer = 0.0;
        denom = 0.0;
        for (t = 0; t < T[n]; t++)
        {

            if (isNaN[n][t] == 0)
            {
//numer = numer + gamma[t][i] * observations[n][t][d]; //indikator
                numer = numer
                                                  //IndicatorFct(observations[n][t][obs_d]);
                    + gamma[t][i] * observations[n][t][obs_d];
                denom = denom + gamma[t][i];
            }
//	printf("obs[%d][%d]=%f, gamma=%f\n", t, obs_d, observations[n][t][obs_d],gamma[t][i]);
        }
        this->updateNumeratorLambda[d] += 1 / Pk[n] * numer;
        this->updateDenominatorLambda[d] += 1 / Pk[n] * denom;

    }
//	printf("\n Poisson type %f/%f\n", this->updateDenominatorLambda[0], this->updateNumeratorLambda[0] );
}


void Poisson::updateAuxiliariesCoupled(double*** observations, double** gamma,
double* Pk, int* T, int n, int i, int statecouple, int** isNaN)
{
    int t, d, l, obs_d;
    double numer, denom;
//printf("BERNOULLI o.revob s=%d, c=%d\n", i, statecouple);
    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        obs_d = OBSERVATIONPOS(d, this->emissionParams->getStart());
//printf("obs_d %d \n", obs_d );
        numer = 0.0;
        denom = 0.0;
        for (t = 0; t < T[n]; t++)
        {
            if (isNaN[n][t] == 0)
            {
                numer = numer
                    + (gamma[t][i] + gamma[t][statecouple])
                    * observations[n][t][obs_d];  //IndicatorFct(observations[n][t][obs_d]);
                denom = denom + (gamma[t][i] + gamma[t][statecouple]);
            }
        }
        this->updateNumeratorLambda[d] += 1 / Pk[n] * numer;
        this->updateDenominatorLambda[d] += 1 / Pk[n] * denom;
    }

}


void Poisson::updateAuxiliariesCoupledRevop(double*** observations,
double** gamma, double* Pk, int* T, int n, int i, int statecouple,
int* state2flag, int* revop, int** isNaN)
{
    int t, d, l, obs_d;
    double numer, denom;

//printf("BERNOULLI s=%d, c=%d\n", i, statecouple);
    for (d = 0; d < this->emissionParams->getD(); d++)
    {

        obs_d = OBSERVATIONPOS(d, this->emissionParams->getStart());
//printf("d= %d  obs_d=%d revObs_d = %d \n ", d, obs_d, revop[obs_d] );
        numer = 0.0;
        denom = 0.0;
        for (t = 0; t < T[n]; t++)
        {
            if (isNaN[n][t] == 0)
            {
                if (state2flag[statecouple] == 1)
                {
                                                  //IndicatorFct(observations[n][t][obs_d])
                    numer = numer + gamma[t][i] *observations[n][t][obs_d]
                        + gamma[t][statecouple]
                                                  //IndicatorFct(
                        * observations[n][t][revop[obs_d]];
//observations[n][t][revop[obs_d]]);
                }
                else
                {
                    numer = numer
                                                  //IndicatorFct(observations[n][t][revop[obs_d]])
                        + gamma[t][i] * observations[n][t][revop[obs_d]]
                        + gamma[t][statecouple]
                        * observations[n][t][obs_d];//IndicatorFct(observations[n][t][obs_d]);
                }
                denom = denom + (gamma[t][i] + gamma[t][statecouple]);
            }
        }
        this->updateNumeratorLambda[d] += 1 / Pk[n] * numer;
        this->updateDenominatorLambda[d] += 1 / Pk[n] * denom;
    }

}


void Poisson::updateCoupledRevop(double ***observations, double* Pk,
int statecouple, int* state2flag, int* revop, double** revGammaAux,
int** isNaN, SEXP emissionPrior, int currN, int ncores)
{

//printf("************ update Coupled REvop\n");
    int skipcounter = 0;
    int d;
    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        this->emissionParams->setPoissonLambda(
            this->updateNumeratorLambda[d] / this->updateDenominatorLambda[d]);

        this->updateNumeratorLambda[d] = 0;
        this->updateDenominatorLambda[d] = 0;

    }
//printf("\n");
}


// was mach ich aus dem prior
double Poisson::Prior(SEXP hyperparams)
{
    return 0;
}


void Poisson::update(double ***observations, double* Pk, int** isNaN,
SEXP emissionPrior, int currN, int ncores)
{
//printf("Entered first update fct from Poisson\n");
    int skipcounter = 0;
    int d;
    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        this->emissionParams->setPoissonLambda(
            this->updateNumeratorLambda[d] / this->updateDenominatorLambda[d]);

        this->updateNumeratorLambda[d] = 0;
        this->updateDenominatorLambda[d] = 0;

    }
//printf("NEW Poisson prob = %f\n", this->emissionParams->getPoissonLambda());

//printf("\n");

}


void Poisson::setParsToTwin(EmissionFunction* myTwinEmission, int currN, double*** observations)
{
    this->emissionParams->setPoissonLambda(myTwinEmission->getParameter()->getPoissonLambda());
    this->updateNumeratorLambda[0] = 0;
    this->updateDenominatorLambda[0] = 0;
}
