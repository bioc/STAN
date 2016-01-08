#include "Bernoulli.h"
#include <cmath>
#include <list>
using namespace std;

#define OBSERVATIONPOS(index, start) (start[index])

Bernoulli::Bernoulli(ParamContainerEmissions *emissionParams)
{

    this->emissionParams = emissionParams;

    int mem = sizeof(double) * this->emissionParams->getD() * 4
        + sizeof(double*) * this->emissionParams->getD() * 2;
    if (DEBUG_MEMORY)
    {
        printf("new->Bernoulli; (%d bytes) \n", mem);
    }

    this->updateNumeratorP = (double*) malloc(
        sizeof(double) * this->emissionParams->getD());
    this->updateDenominatorP = (double*) malloc(
        sizeof(double) * this->emissionParams->getD());

    int d1, d2;
    for (d1 = 0; d1 < this->emissionParams->getD(); d1++)
    {
        this->updateNumeratorP[d1] = 0;
        this->updateDenominatorP[d1] = 0;
    }
}


Bernoulli::~Bernoulli()
{

    int d1;
    free(this->updateNumeratorP);
    free(this->updateDenominatorP);

    int mem = sizeof(double) * this->emissionParams->getD() * 4
        + sizeof(double*) * this->emissionParams->getD() * 2;
    if (DEBUG_MEMORY)
    {
        printf("delete->Bernoulli; (%d bytes) ", mem);
    }
//delete this->emissionParams;
}


double Bernoulli::getP()
{
    return this->emissionParams->getBernoulliP();
}


double Bernoulli::calcEmissionProbability(double *obs, int isna, int currN)
{

    int i, j, obs_i;
    double probability = 1;
    if (isna == 0)
    {
        for (i = 0; i < this->emissionParams->getD(); i++)
        {
//	printf("start=%d, D=%d\n",this->emissionParams->getD(), this->emissionParams->getStart()[i]);

            obs_i = OBSERVATIONPOS(i, this->emissionParams->getStart());

            if (obs[obs_i] != obs[obs_i])         //check NAN
            {
                isna = 1;
                break;
            }
            else
            {

                probability = probability
                    * pow(this->emissionParams->getBernoulliP(), obs[obs_i])
                    * pow((1 - this->emissionParams->getBernoulliP()),
                    (1 - obs[obs_i]));

            }
        }

    }

    if (probability < 1e-100)
    {
        probability = 1e-100;
    }
    return probability;
}


void Bernoulli::updateAuxiliaries(double*** observations, double** gamma,
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
        this->updateNumeratorP[d] += 1 / Pk[n] * numer;
        this->updateDenominatorP[d] += 1 / Pk[n] * denom;

    }
//	printf("\n Bernoulli type %f/%f\n", this->updateDenominatorP[0], this->updateNumeratorP[0] );
}


int Bernoulli::IndicatorFct(int obs)
{
    if (obs == 1)
        return 1;
    else
        return 0;
}


void Bernoulli::updateAuxiliariesCoupled(double*** observations, double** gamma,
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
        this->updateNumeratorP[d] += 1 / Pk[n] * numer;
        this->updateDenominatorP[d] += 1 / Pk[n] * denom;
    }

}


void Bernoulli::updateAuxiliariesCoupledRevop(double*** observations,
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
        this->updateNumeratorP[d] += 1 / Pk[n] * numer;
        this->updateDenominatorP[d] += 1 / Pk[n] * denom;
    }

}


void Bernoulli::updateCoupledRevop(double ***observations, double* Pk,
int statecouple, int* state2flag, int* revop, double** revGammaAux,
int** isNaN, SEXP emissionPrior, int currN, int ncores)
{

//printf("************ update Coupled REvop\n");
    int skipcounter = 0;
    int d;
    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        this->emissionParams->setBernoulliP(
            this->updateNumeratorP[d] / this->updateDenominatorP[d]);

        this->updateNumeratorP[d] = 0;
        this->updateDenominatorP[d] = 0;

    }
//printf("\n");
}


// was mach ich aus dem prior
double Bernoulli::Prior(SEXP hyperparams)
{
/*int i, j;
for (i = 0; i < this->emissionParams->getD(); i++) {
    for (j = 0; j < this->emissionParams->getD(); j++) {
        REAL(getListElement(hyperparams, "cov"))[i
                + j * this->emissionParams->getD()] =
                this->emissionParams->getGaussianSIGMA()[i][j];
        // Rprintf("%f ", this->emissionParams->getGaussianSIGMA()[i][j]);
    }
    //Rprintf("\n");
}

// call solnp from R for optimization
SEXP call = PROTECT(lang2(install("calldiwish"), hyperparams));
SEXP res = PROTECT(eval(call, R_GlobalEnv));
double out = REAL(res)[0];
//Rprintf("%f\n", out);
UNPROTECT(2);
//return out;*/
    return 0;
}


void Bernoulli::update(double ***observations, double* Pk, int** isNaN,
SEXP emissionPrior, int currN, int ncores)
{
//printf("Entered first update fct from Bernoulli\n");
    int skipcounter = 0;
    int d;
    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        this->emissionParams->setBernoulliP(
            this->updateNumeratorP[d] / this->updateDenominatorP[d]);

        this->updateNumeratorP[d] = 0;
        this->updateDenominatorP[d] = 0;

    }
//printf("NEW Bernoulli prob = %f\n", this->emissionParams->getBernoulliP());

//printf("\n");

}


void Bernoulli::setParsToTwin(EmissionFunction* myTwinEmission, int currN, double*** observations)
{
    this->emissionParams->setBernoulliP(myTwinEmission->getParameter()->getBernoulliP());
    int d;
    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        this->updateNumeratorP[d] = 0;
        this->updateDenominatorP[d] = 0;

    }
}
