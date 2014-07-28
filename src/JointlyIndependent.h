
#ifndef JointlyIndependent_HEADER
#define JointlyIndependent_HEADER

#include "EmissionFunction.h"
#include "ParamContainerEmissions.h"
#include "DebugConstants.h"
#include "RAccessUtils.h"

#include <list>
using namespace std;

class JointlyIndependent : public EmissionFunction
{
protected:

    list<EmissionFunction*> listEmissions;

public:
    JointlyIndependent(list<EmissionFunction*> efb, ParamContainerEmissions *emissionParams);
    JointlyIndependent() {}
    ~JointlyIndependent();
    double calcEmissionProbability(double *obs, int isNaN);
    virtual void updateAuxiliaries(double*** observations, double** gamma, double* Pk, int* T, int n, int i, int** isNaN);
    virtual void updateAuxiliariesCoupled(double*** observations, double** gamma, double* Pk, int* T, int n, int i, int statecouple, int** isNaN);
    virtual void updateAuxiliariesCoupledRevop(double*** observations, double** gamma, double* Pk, int* T, int n, int i, int statecouple, int* state2flag, int* revop, int** isNaN);
    virtual void updateCoupledRevop(double ***observations, double* Pk, int statecouple, int* state2flag, int* revop, double** revGammaAux, int** isNaN, SEXP emissionPrior, int currN);
    virtual void update(double ***observations, double* Pk, int** isNaN, SEXP emissionPrior, int currN);
    std::list<EmissionFunction*> getEmissionFunctionList();
    double Prior(SEXP hyperparams);

private:
    int IndicatorFct(int val);
};
#endif
