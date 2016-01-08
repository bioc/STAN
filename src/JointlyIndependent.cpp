#include "JointlyIndependent.h"
//#include <iostream>
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
    std::list<EmissionFunction*>::iterator pos;
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

    std::list<EmissionFunction*>::iterator pos;
    for (pos = this->listEmissions.begin(); pos != this->listEmissions.end();
        pos++)
    {
        mem = mem + sizeof(double) * (*pos)->getParameter()->getD() * 3;
        EmissionFunction* currEmission = *pos;
//	printf("%d\n", (currEmission)->getParameter()->getWhichOne());
//	 listEmissions.erase(pos);
        delete currEmission;
//listEmissions.remove(*pos);
//	pos++;
    }
/*
    printf("D=%d\n", D);
    int d;
    for(d=0; d<2; d++) {
        printf("%d\n", d);
        printf("%d\n", emissionArray[d]->getParameter()->getWhichOne());
        delete emissionArray[d];
    }
*/
    this->listEmissions.clear();
    this->ContainerList.clear();
    delete this->emissionParams;
    if (DEBUG_MEMORY)
    {
        printf("delete->JointlyIndependent; (%d bytes) ", mem);
    }

}


double JointlyIndependent::calcEmissionProbability(double *obs, int isna, int currN)
{
    double probalility = 1;
    std::list<EmissionFunction*>::iterator pos;
    for (pos = this->listEmissions.begin(); pos != this->listEmissions.end();
        pos++)
    {

        probalility = probalility * (*pos)->calcEmissionProbability(obs, isna, currN);
    }
    if(probalility > 1e30)
    {
//Rprintf("ji-prob=%.100f\n", probalility);

    }
    if (probalility < 1e-300)
    {
//Rprintf("ji=>small\n");

        probalility = 1e-300;
    }

    return probalility;
}


void JointlyIndependent::updateAuxiliaries(double*** observations,
double** gamma, double* Pk, int* T, int n, int i, int** isNaN)
{
//	printf("Beginn updateAuxiliaries %d\n", n);
    int t;
    for(t=0; t<T[n]; t++)
    {
        this->emissionParams->setGammaAux(gamma[t][i], n, t);
    }

    std::list<EmissionFunction*>::iterator pos;
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
    int t;
    for(t=0; t<T[n]; t++)
    {
        this->emissionParams->setGammaAux(gamma[t][i], n, t);
    }
    std::list<EmissionFunction*>::iterator pos;
    for (pos = this->listEmissions.begin(); pos != this->listEmissions.end();
        pos++)
    {
//printf("UPDATE WHICH: %d \n", (*pos)->getParameter()->getWhichOne());
        (*pos)->updateAuxiliariesCoupled(observations, gamma, Pk, T, n, i,
            statecouple, isNaN);
    }

//printf("\n");
}


void JointlyIndependent::updateAuxiliariesCoupledRevop(double*** observations,
double** gamma, double* Pk, int* T, int n, int i, int statecouple,
int* state2flag, int* revop, int** isNaN)
{
    std::list<EmissionFunction*>::iterator pos;

    int t;
    for(t=0; t<T[n]; t++)
    {
        this->emissionParams->setGammaAux(gamma[t][i], n, t);
    }
    for (pos = this->listEmissions.begin(); pos != this->listEmissions.end();
        pos++)
    {
//	printf("UPDATE WHICH: %d \n", (*pos)->getParameter()->getWhichOne());
        (*pos)->updateAuxiliariesCoupledRevop(observations, gamma, Pk, T, n, i,
            statecouple, state2flag, revop, isNaN);

    }

//printf("\n");

}


void JointlyIndependent::updateCoupledRevop(double ***observations, double* Pk,
int statecouple, int* state2flag, int* revop,
double** revGammaAux, int** isNaN, SEXP emissionPrior, int currN, int ncores)
{

    std::list<EmissionFunction*>::iterator pos;

    for (pos = this->listEmissions.begin(); pos != this->listEmissions.end();
        pos++)
    {

        (*pos)->updateCoupledRevop(observations, Pk, statecouple, state2flag, revop, revGammaAux, isNaN, emissionPrior, currN, ncores);

    }

}


// was mach ich aus dem prior
double JointlyIndependent::Prior(SEXP hyperparams)
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
    SEXP call = PROTECT(lang2(install("calldiwish"), hyperparams));
    SEXP res = PROTECT(eval(call, R_GlobalEnv));
    double out = REAL(res)[0];
//Rprintf("%f\n", out);
    UNPROTECT(2);
    return out;

}


void JointlyIndependent::update(double ***observations, double* Pk, int** isNaN,
SEXP emissionPrior, int currN, int ncores)
{

    std::list<EmissionFunction*>::iterator pos;
//int len = listEmissions.size();
//EmissionFunction **alldims = (EmissionFunction**)malloc(sizeof(EmissionFunction*)*len);
//int i = 0;
    for (pos = this->listEmissions.begin(); pos != this->listEmissions.end();
        pos++)
    {
//printf("UPDATE WHICH: %d \n", (*pos)->getParameter()->getWhichOne());
        (*pos)->update(observations, Pk, isNaN, emissionPrior, currN, ncores);
//alldims[i] = (*pos);
//i++;
    }

//#pragma omp parallel for
//for(i=0; i<len; i++) {
//	alldims[i]->update(observations, Pk, isNaN, emissionPrior, currN, ncores);
//}

}


void JointlyIndependent::setParsToTwin(EmissionFunction* myTwinEmission, int currN, double*** observations)
{

    std::list<EmissionFunction*>::iterator pos;
    std::list<EmissionFunction*>::iterator posTwin;
    int i=0;
//posTwin++;
//EmissionFunction* currEmission = (*posTwin);
    std::list<EmissionFunction*> myTwinEmissionList =  myTwinEmission->getEmissionFunctionList();
    int currSize = myTwinEmissionList.size();
    EmissionFunction** myemissions = (EmissionFunction**)malloc(sizeof(EmissionFunction*)*currSize);
    for (posTwin = myTwinEmissionList.begin(); posTwin != myTwinEmissionList.end();
        posTwin++)
    {
        myemissions[i] = (*posTwin);
//Rprintf("type=%d\n", myemissions[i]->getParameter()->getWhichOne());
        i++;
    }
//printf("i=%d\n", currSize);
//myemissions[0] = myTwinEmissionList.front();
    i=0;
    for (pos = this->listEmissions.begin(); pos != this->listEmissions.end();
        pos++)
    {
        (*pos)->setParsToTwin(myemissions[i], currN, observations);
        i++;
    }
    free(myemissions);
}


list<EmissionFunction*> JointlyIndependent::getEmissionFunctionList()
{
    return this->listEmissions;
}


void JointlyIndependent::computeShared(EmissionFunction** myEmissions, int nStates)
{
    int k,d;
    int nEmissions = this->listEmissions.size();
    EmissionFunction*** currEmissions = (EmissionFunction***)malloc(sizeof(EmissionFunction**)*nEmissions);

    for(d=0; d<nEmissions; d++)
    {
        currEmissions[d] = (EmissionFunction**)malloc(sizeof(EmissionFunction*)*nStates);
    }

    std::list<EmissionFunction*> listEF;
    for (k = 0; k < nStates; k++)
    {
        std::list<EmissionFunction*>::iterator pos;
        listEF = myEmissions[k]->getEmissionFunctionList();
        d = 0;
        for (pos = listEF.begin(); pos!= listEF.end(); pos++)
        {
            currEmissions[d][k] = (*pos);
//printf("wo=%d\n", (*pos)->getParameter()->getWhichOne());
            d++;
        }
    }

/*
        d=0;
        list<EmissionFunction*>::iterator pos;
        for (pos = this->listEmissions.begin(); pos != this->listEmissions.end(); pos++) {
            for(k=0; k<nStates; k++) {
                currEmissions[d][k] = (*pos);
            }
            d++;
        }*/
    d=0;
    std::list<EmissionFunction*>::iterator pos;
    for (pos = this->listEmissions.begin(); pos != this->listEmissions.end(); pos++)
    {
        (*pos)->computeShared(currEmissions[d], nStates);
        d++;
    }

/*for(d=0; d<nEmissions; d++) {
    this->computeShared(currEmissions, nStates);
}*/
}


void JointlyIndependent::resetShared()
{
    list<EmissionFunction*>::iterator pos;

    for (pos = this->listEmissions.begin(); pos != this->listEmissions.end();
        pos++)
    {

        (*pos)->resetShared();

    }

}
