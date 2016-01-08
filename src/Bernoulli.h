
#ifndef BERNOULLI_HEADER
#define BERNOULLI_HEADER

#include "EmissionFunction.h"
#include "ParamContainerEmissions.h"
#include "DebugConstants.h"
#include "RAccessUtils.h"
#include <list>
using namespace std;

class Bernoulli : public EmissionFunction
{
    protected:
        double* updateNumeratorP;
        double* updateDenominatorP;
        double* updateNumeratorMU;
        double* updateDenominatorMU;
        double*** gammaAux;
        double** updateNumeratorSIGMA;
        double** updateDenominatorSIGMA;

    public:
        Bernoulli(ParamContainerEmissions *emissionParams);
        Bernoulli() {}
        ~Bernoulli();
        double calcEmissionProbability(double *obs, int isNaN, int currN);
        double getP();
        virtual void updateAuxiliaries(double*** observations, double** gamma, double* Pk, int* T, int n, int i, int** isNaN);
        virtual void updateAuxiliariesCoupled(double*** observations, double** gamma, double* Pk, int* T, int n, int i, int statecouple, int** isNaN);
        virtual void updateAuxiliariesCoupledRevop(double*** observations, double** gamma, double* Pk, int* T, int n, int i, int statecouple, int* state2flag, int* revop, int** isNaN);
        virtual void updateCoupledRevop(double ***observations, double* Pk, int statecouple, int* state2flag, int* revop, double** revGammaAux, int** isNaN, SEXP emissionPrior, int currN, int ncores);
        virtual void update(double ***observations, double* Pk, int** isNaN, SEXP emissionPrior, int currN, int ncores);
        double Prior(SEXP hyperparams);
        virtual void setParsToTwin(EmissionFunction* myTwinEmission, int currN, double*** observations);
    private:
        int IndicatorFct(int val);
};
#endif
