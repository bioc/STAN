#include "JointlyIndependent.h"
#include <cmath>
#include <list>
using namespace std;

#define ARRAYSIZE(a) (sizeof(a) / sizeof(a[0]))
#define COLSIZE(a) ARRAYSIZE(a[0])

JointlyIndependent::JointlyIndependent(list<EmissionFunction*> efb, ParamContainerEmissions *emissionParams)
{

    this->emissionParams = emissionParams;
    this->listEmissions = efb;
    this->ContainerList = efb;

    int mem = 0;
    list<EmissionFunction*>::iterator pos;
    for (pos = this->listEmissions.begin(); pos != this->listEmissions.end();
            pos++)
    {
        mem = mem + sizeof(double) * (*pos)->getParameter()->getD() * 3;

    }

    if (DEBUG_MEMORY)
    {
        printf("new->JointlyIndependent; (%d bytes) \n", mem);
    }

}


JointlyIndependent::~JointlyIndependent()
{
    int mem = 0;

    list<EmissionFunction*>::iterator pos;
    for (pos = this->listEmissions.begin(); pos != this->listEmissions.end();
            pos++)
    {
        mem = mem + sizeof(double) * (*pos)->getParameter()->getD() * 3;

    }

    if (DEBUG_MEMORY)
    {
        printf("delete->JointlyIndependent; (%d bytes) ", mem);
    }

}


double JointlyIndependent::calcEmissionProbability(double *obs, int isna)
{
    double probalility = 1;
    list<EmissionFunction*>::iterator pos;
    for (pos = this->listEmissions.begin(); pos != this->listEmissions.end();
            pos++)
    {
        probalility = probalility * (*pos)->calcEmissionProbability(obs, isna);

    }
    if (probalility < 1e-100)
    {
        probalility = 1e-100;
    }

    return probalility;
}


void JointlyIndependent::updateAuxiliaries(double*** observations,
        double** gamma, double* Pk, int* T, int n, int i, int** isNaN)
{
    list<EmissionFunction*>::iterator pos;
    for (pos = this->listEmissions.begin(); pos != this->listEmissions.end();
            pos++)
    {

        (*pos)->updateAuxiliaries(observations, gamma, Pk, T, n, i, isNaN);
    }

}


int JointlyIndependent::IndicatorFct(int obs)
{
    if (obs == 1)
        return 1;
    else
        return 0;
}


void JointlyIndependent::updateAuxiliariesCoupled(double*** observations,
        double** gamma, double* Pk, int* T, int n, int i, int statecouple,
        int** isNaN)
{
    list<EmissionFunction*>::iterator pos;
    for (pos = this->listEmissions.begin(); pos != this->listEmissions.end();
            pos++)
    {
        (*pos)->updateAuxiliariesCoupled(observations, gamma, Pk, T, n, i,
                                         statecouple, isNaN);
    }
    int t;
}


void JointlyIndependent::updateAuxiliariesCoupledRevop(double*** observations,
        double** gamma, double* Pk, int* T, int n, int i, int statecouple,
        int* state2flag, int* revop, int** isNaN)
{
    list<EmissionFunction*>::iterator pos;

    for (pos = this->listEmissions.begin(); pos != this->listEmissions.end();
            pos++)
    {
        (*pos)->updateAuxiliariesCoupledRevop(observations, gamma, Pk, T, n, i,
                                              statecouple, state2flag, revop, isNaN);

    }
    int t;
}


void JointlyIndependent::updateCoupledRevop(double ***observations, double* Pk,
        int statecouple, int* state2flag, int* revop,
        double** revGammaAux, int** isNaN, SEXP emissionPrior, int currN)
{

    list<EmissionFunction*>::iterator pos;

    for (pos = this->listEmissions.begin(); pos != this->listEmissions.end();
            pos++)
    {

        (*pos)->updateCoupledRevop(observations, Pk, statecouple, state2flag, revop, revGammaAux, isNaN, emissionPrior, currN);

    }

}


double JointlyIndependent::Prior(SEXP hyperparams)
{
    double out = 1;
    return out;
}


void JointlyIndependent::update(double ***observations, double* Pk, int** isNaN,
                                SEXP emissionPrior, int currN)
{

    list<EmissionFunction*>::iterator pos;
    for (pos = this->listEmissions.begin(); pos != this->listEmissions.end();
            pos++)
    {
        (*pos)->update(observations, Pk, isNaN, emissionPrior, currN);

    }

}


list<EmissionFunction*> JointlyIndependent::getEmissionFunctionList()
{
    return this->listEmissions;
}
