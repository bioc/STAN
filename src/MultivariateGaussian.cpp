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


double MultivariateGaussian::calcEmissionProbability(double *obs, int isna)
{
    int i, j;
    double exponent = 0;
    double first_term = 0;
    int obs_i = 0;
    int obs_j = 0;
    first_term = pow(2.5066282746310002, this->emissionParams->getD())
                 * sqrt(this->emissionParams->getGaussianDET());
    double density = 1;

    if (isna == 0)
    {
        for (i = 0; i < this->emissionParams->getD(); i++)
        {
            obs_i = OBSERVATIONPOS(i, this->emissionParams->getStart());
            if (obs[obs_i] != obs[obs_i])
            {
                isna = 1;
                break;
            }
            for (j = 0; j < this->emissionParams->getD(); j++)
            {
                obs_j = OBSERVATIONPOS(j, this->emissionParams->getStart());
                if (obs[obs_j] != obs[obs_j])
                {
                    isna = 1;
                    break;
                }
                exponent += (obs[obs_i]
                             - this->emissionParams->getGaussianMU()[i][0])
                            * this->emissionParams->getGaussianINVSIGMA()[i][j]
                            * (obs[obs_j]
                               - this->emissionParams->getGaussianMU()[j][0]);
            }
        }
        density = exp(-0.5 * exponent) / first_term;
    }
    if (density > 1e20)
    {
        error("Ill-conditioned covariance matrix!\n");
    }
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

    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        numer = 0.0;
        denom = 0.0;
        obs_d = OBSERVATIONPOS(d, this->emissionParams->getStart());
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
        double** revGammaAux, int** isNaN, SEXP emissionPrior, int currN)
{

    int skipcounter = 0;
    int d;
    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        this->emissionParams->setGaussianMUelement(
            this->updateNumeratorMU[d] / this->updateDenominatorMU[d], d);
        this->updateNumeratorMU[d] = 0;
        this->updateDenominatorMU[d] = 0;
    }

    double numer, denom;

    int d1, d2, t, l, obs_d1, obs_d2;
    int n;

    int lower_n = 0;
    int upper_n = this->emissionParams->getNsample();
    if(currN != -1)
    {
        lower_n = currN;
        upper_n = currN+1;
    }

    for (n = lower_n; n < upper_n; n++)
    {

        for (d1 = 0; d1 < this->emissionParams->getD(); d1++)
        {

            obs_d1 = OBSERVATIONPOS(d1, this->emissionParams->getStart());
            for (d2 = d1; d2 < this->emissionParams->getD(); d2++)
            {
                obs_d2 = OBSERVATIONPOS(d2, this->emissionParams->getStart());
                numer = 0.0;
                denom = 0.0;
                for (t = 0; t < this->emissionParams->getT()[n]; t++)
                {

                    if (isNaN[n][t] == 0)
                    {

                        if (this->emissionParams->getGammaAux()[n][t] < 1e-5)
                        {
                            skipcounter++;
                        }

                        if (state2flag[statecouple] == 1)
                        {
                            numer =
                                numer
                                + this->emissionParams->getGammaAux()[n][t]
                                * (observations[n][t][obs_d1]
                                   - this->emissionParams->getGaussianMU()[d1][0])
                                * (observations[n][t][obs_d2]
                                   - this->emissionParams->getGaussianMU()[d2][0])
                                + revGammaAux[n][t]
                                * (observations[n][t][revop[obs_d1]]
                                   - this->emissionParams->getGaussianMU()[d1][0])
                                * (observations[n][t][revop[obs_d2]]
                                   - this->emissionParams->getGaussianMU()[d2][0]);
                        }
                        else
                        {
                            numer =
                                numer
                                + this->emissionParams->getGammaAux()[n][t]
                                * (observations[n][t][revop[obs_d1]]
                                   - this->emissionParams->getGaussianMU()[d1][0])
                                * (observations[n][t][revop[obs_d2]]
                                   - this->emissionParams->getGaussianMU()[d2][0])
                                + revGammaAux[n][t]
                                * (observations[n][t][obs_d1]
                                   - this->emissionParams->getGaussianMU()[d1][0])
                                * (observations[n][t][obs_d2]
                                   - this->emissionParams->getGaussianMU()[d2][0]);
                        }

                        denom = denom
                                + this->emissionParams->getGammaAux()[n][t]
                                + revGammaAux[n][t];
                    }

                }

                if (LENGTH(emissionPrior) > 0)
                {
                    this->updateNumeratorSIGMA[d1][d2] += REAL(
                            coerceVector(getListElement(emissionPrior, "S"),
                                         REALSXP))[d1
                                                          + d2 * this->emissionParams->getD()];
                    this->updateDenominatorSIGMA[d1][d2] += REAL(
                            getListElement(emissionPrior, "v"))[0]
                                                            + this->emissionParams->getD() + 1;
                }

                this->updateNumeratorSIGMA[d1][d2] += 1.0 / Pk[n] * numer;
                this->updateDenominatorSIGMA[d1][d2] += 1.0 / Pk[n] * denom;
                if (d1 != d2)
                {
                    this->updateNumeratorSIGMA[d2][d1] =
                        this->updateNumeratorSIGMA[d1][d2];
                    this->updateDenominatorSIGMA[d2][d1] =
                        this->updateDenominatorSIGMA[d1][d2];
                }
            }

        }
    }

    for (d1 = 0; d1 < this->emissionParams->getD(); d1++)
    {
        for (d2 = 0; d2 < this->emissionParams->getD(); d2++)
        {
            this->emissionParams->setGaussianSIGMAelement(
                this->updateNumeratorSIGMA[d1][d2]
                / this->updateDenominatorSIGMA[d1][d2], d1, d2);
            this->emissionParams->setGaussianINVSIGMAelement(
                this->updateNumeratorSIGMA[d1][d2]
                / this->updateDenominatorSIGMA[d1][d2], d1, d2);
            this->updateNumeratorSIGMA[d1][d2] = 0;
            this->updateDenominatorSIGMA[d1][d2] = 0;
        }
    }

    inverse(this->emissionParams->getGaussianINVSIGMA(),
            this->emissionParams->getD());
    this->emissionParams->setGaussianDET(
        matrixDet(this->emissionParams->getGaussianSIGMA(),
                  this->emissionParams->getD()));

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
        }
    }

    // call solnp from R for optimization
    SEXP call = PROTECT(lang2( getListElement(hyperparams, "calldiwish"), hyperparams));
    SEXP res = PROTECT(eval(call, R_GlobalEnv));
    double out = REAL(res)[0];
    UNPROTECT(2);
    return out;

}


void MultivariateGaussian::update(double ***observations, double* Pk,
                                  int** isNaN, SEXP emissionPrior, int currN)
{
    int skipcounter = 0;
    int d;
    for (d = 0; d < this->emissionParams->getD(); d++)
    {
        this->emissionParams->setGaussianMUelement(
            this->updateNumeratorMU[d] / this->updateDenominatorMU[d], d);
        this->updateNumeratorMU[d] = 0;
        this->updateDenominatorMU[d] = 0;
    }

    double numer, denom;
    int d1, d2, t, n, l, obs_d1, obs_d2;
    int lower_n = 0;
    int upper_n = this->emissionParams->getNsample();
    if(currN != -1)
    {
        lower_n = currN;
        upper_n = currN+1;
    }

    for (n = lower_n; n < upper_n; n++)
    {
        for (d1 = 0; d1 < this->emissionParams->getD(); d1++)
        {
            obs_d1 = OBSERVATIONPOS(d1, this->emissionParams->getStart());
            for (d2 = d1; d2 < this->emissionParams->getD(); d2++)
            {
                obs_d2 = OBSERVATIONPOS(d2, this->emissionParams->getStart());

                numer = 0.0;
                denom = 0.0;
                for (t = 0; t < this->emissionParams->getT()[n]; t++)
                {
                    if (isNaN[n][t] == 0)
                    {
                        if (this->emissionParams->getGammaAux()[n][t] < 1e-5)
                        {
                            skipcounter++;
                        }

                        numer =
                            numer
                            + this->emissionParams->getGammaAux()[n][t]
                            * (observations[n][t][obs_d1]
                               - this->emissionParams->getGaussianMU()[d1][0])
                            * (observations[n][t][obs_d2]
                               - this->emissionParams->getGaussianMU()[d2][0]);
                        denom = denom
                                + this->emissionParams->getGammaAux()[n][t];
                    }
                }

                if (LENGTH(emissionPrior) > 0)
                {
                    this->updateNumeratorSIGMA[d1][d2] += REAL(
                            coerceVector(getListElement(emissionPrior, "S"),
                                         REALSXP))[d1
                                                          + d2 * this->emissionParams->getD()];
                    this->updateDenominatorSIGMA[d1][d2] += REAL(
                            getListElement(emissionPrior, "v"))[0]
                                                            + this->emissionParams->getD() + 1;
                }

                this->updateNumeratorSIGMA[d1][d2] += 1.0 / Pk[n] * numer;
                this->updateDenominatorSIGMA[d1][d2] += 1.0 / Pk[n] * denom;
                if (d1 != d2)
                {
                    this->updateNumeratorSIGMA[d2][d1] =
                        this->updateNumeratorSIGMA[d1][d2];
                    this->updateDenominatorSIGMA[d2][d1] =
                        this->updateDenominatorSIGMA[d1][d2];
                }
            }
        }
    }

    for (d1 = 0; d1 < this->emissionParams->getD(); d1++)
    {
        for (d2 = 0; d2 < this->emissionParams->getD(); d2++)
        {
            this->emissionParams->setGaussianSIGMAelement(
                this->updateNumeratorSIGMA[d1][d2]
                / this->updateDenominatorSIGMA[d1][d2], d1, d2);
            this->emissionParams->setGaussianINVSIGMAelement(
                this->updateNumeratorSIGMA[d1][d2]
                / this->updateDenominatorSIGMA[d1][d2], d1, d2);
            this->updateNumeratorSIGMA[d1][d2] = 0;
            this->updateDenominatorSIGMA[d1][d2] = 0;
        }
    }
    inverse(this->emissionParams->getGaussianINVSIGMA(),
            this->emissionParams->getD());
    this->emissionParams->setGaussianDET(
        matrixDet(this->emissionParams->getGaussianSIGMA(),
                  this->emissionParams->getD()));

}
