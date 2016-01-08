
#ifndef MULTGAUSS_HEADER
#define MULTGAUSS_HEADER

#include "EmissionFunction.h"
#include "ParamContainerEmissions.h"
#include "DebugConstants.h"
#include "RAccessUtils.h"

class MultivariateGaussian : public EmissionFunction
{
    protected:
        double* updateNumeratorMU;
        double* updateDenominatorMU;
        double*** gammaAux;

    public:
        MultivariateGaussian(ParamContainerEmissions *emissionParams);
        MultivariateGaussian() {}
        ~MultivariateGaussian();
        double calcEmissionProbability(double *obs, int isNaN, int currN);
        virtual void updateAuxiliaries(double*** observations, double** gamma, double* Pk, int* T, int n, int i, int** isNaN);
        virtual void updateAuxiliariesCoupled(double*** observations, double** gamma, double* Pk, int* T, int n, int i, int statecouple, int** isNaN);
        virtual void updateAuxiliariesCoupledRevop(double*** observations, double** gamma, double* Pk, int* T, int n, int i, int statecouple, int* state2flag, int* revop, int** isNaN);
        virtual void updateCoupledRevop(double ***observations, double* Pk, int statecouple, int* state2flag, int* revop, double** revGammaAux, int** isNaN, SEXP emissionPrior, int currN, int ncores);
        virtual void update(double ***observations, double* Pk, int** isNaN, SEXP emissionPrior, int currN, int ncores);
        virtual double** getUpdateNumeratorSigma();
        virtual double** getUpdateDenominatorSigma();
        virtual void computeShared(EmissionFunction** myEmissions, int nStates);
        virtual void resetShared();
        double Prior(SEXP hyperparams);
        virtual void setParsToTwin(EmissionFunction* myTwinEmission, int currN, double*** observations);
};
#endif
