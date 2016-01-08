#include "NegativeBinomial.h"
#include <cmath>
#define OBSERVATIONPOS(index, start) (start[index])

NegativeBinomial::NegativeBinomial(
ParamContainerEmissions *emissionParams)
{

    this->emissionParams = emissionParams;

    int mem = sizeof(double) * this->emissionParams->getD() * 4
        + sizeof(double*) * this->emissionParams->getD() * 2;

    if (DEBUG_MEMORY)
    {
       Rprintf("new->NegativeBinomial; (%d bytes) \n", mem);
    }

/*	this->updateNumeratorMU = (double*) malloc(
                sizeof(double) * this->emissionParams->getD());
    this->updateDenominatorMU = (double*) malloc(
                sizeof(double) * this->emissionParams->getD());

        int d1, d2;
        for (d1 = 0; d1 < this->emissionParams->getD(); d1++) {
            this->updateNumeratorMU[d1] = 0;
            this->updateDenominatorMU[d1] = 0;
        }
*/
}


NegativeBinomial::~NegativeBinomial()
{

    int d1;
    if (DEBUG_MEMORY)
    {
        Rprintf("delete->NegativeBinomial; (%d bytes) ", 0);
    }
    delete this->emissionParams;
}


double NegativeBinomial::calcEmissionProbability(double *obs, int isna, int currN)
{

    int i, obs_i, D;
    double probability = 1;
    int* myStart = this->emissionParams->getStart();
    double myPi = this->emissionParams->getPiNB();
//	printf("sF=%f\n",this->getParameter()->getSizeFactorNB()[currN]);

    D = this->emissionParams->getD();
    if(isna == -1 | isna == 1)
    {
        for (i = 0; i < D; i++)
        {
            obs_i = 0;                            //myStart[i];
            if (obs[obs_i] != obs[obs_i])
            {
                isna = 1;
                break;
            }
//probability = dnbinom_mu(obs[obs_i], this->emissionParams->getSizeNB(), this->emissionParams->getMuNB()/(this->getParameter()->getSizeFactorNB()[currN]), 0);
            if(obs[obs_i] == 0)
            {
                probability = myPi+(1-myPi)*dnbinom_mu(obs[obs_i], this->emissionParams->getSizeNB(), this->emissionParams->getMuNB()/(this->getParameter()->getSizeFactorNB()[currN]), 0);
            }
            else
            {
                probability = (1-myPi)*dnbinom_mu(obs[obs_i], this->emissionParams->getSizeNB(), this->emissionParams->getMuNB()/(this->getParameter()->getSizeFactorNB()[currN]), 0);
            }
        }
    }
    else
    {
        obs_i = myStart[0];
        int myIndex = (int)obs[obs_i];
        probability = this->getParameter()->getUniqueObsProb()[currN][myIndex];
//	probability = (1-myPi)*dnbinom_mu(obs[obs_i], this->emissionParams->getSizeNB(), this->emissionParams->getMuNB()/(this->getParameter()->getSizeFactorNB()[currN]), 0);

    }

    if(probability < 0)
    {
        Rprintf("%f\n", probability);
        error("Negative probability in NegativeBinomial!");

    }

    if (probability < 1e-100)
    {
        probability = 1e-100;
    }

    return probability;
}


void NegativeBinomial::updateAuxiliaries(double*** observations,
double** gamma, double* Pk, int* T, int n, int i, int** isNaN)
{
// TODO !!!
/*int t;
for (t = 0; t < T[n]; t++) {
    this->emissionParams->setGammaAux(gamma[t][i], n, t);
}*/
}


void NegativeBinomial::updateAuxiliariesCoupled(double*** observations,
double** gamma, double* Pk, int* T, int n, int i, int statecouple,
int** isNaN)
{

/*	int t;
        for (t = 0; t < T[n]; t++) {
            this->emissionParams->setGammaAux(gamma[t][i], n, t);
        }*/
}


void NegativeBinomial::updateAuxiliariesCoupledRevop(double*** observations,
double** gamma, double* Pk, int* T, int n, int i, int statecouple,
int* state2flag, int* revop, int** isNaN)
{

/*int t, d, l, obs_d;
double numer, denom;

//printf(" GAUSSIAN s=%d, c=%d\n", i, statecouple);
for (d = 0; d < this->emissionParams->getD(); d++) {
    numer = 0.0;
    denom = 0.0;
    obs_d = OBSERVATIONPOS(d, this->emissionParams->getStart());
    //printf("GAUSSIAN d= %d obs_d=%d \n",d, obs_d);
    for (t = 0; t < T[n]; t++) {
        if (isNaN[n][t] == 0) {

if (state2flag[statecouple] == 1) {
numer = numer + gamma[t][i] * observations[n][t][obs_d]
+ gamma[t][statecouple]
* observations[n][t][revop[obs_d]]; //// d
} else {
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

for (t = 0; t < T[n]; t++) {
this->emissionParams->setGammaAux(gamma[t][i], n, t);
}
*/
/*int t;
    for (t = 0; t < T[n]; t++) {
        this->emissionParams->setGammaAux(gamma[t][i], n, t);
    }*/
}


void NegativeBinomial::updateCoupledRevop(double ***observations,
double* Pk, int statecouple, int* state2flag, int* revop,
double** revGammaAux, int** isNaN, SEXP emissionPrior, int currN, int ncores)
{
    int n;
    SEXP myObs, pars, myGammas, myd;

    int d;
/*		for (d = 0; d < this->emissionParams->getD(); d++) {
            printf("%f ", this->updateNumeratorMU[d] / this->updateDenominatorMU[d]);
            this->updateNumeratorMU[d] = 0;
            this->updateDenominatorMU[d] = 0;

        }
        printf("\n");*/
//Rprintf("Getting ready to optimize!\n");
    PROTECT(pars=NEW_NUMERIC(3));
    NUMERIC_POINTER(pars)[0] = this->emissionParams->getMuNB();
    NUMERIC_POINTER(pars)[1] = this->emissionParams->getSizeNB();
    NUMERIC_POINTER(pars)[2] = this->emissionParams->getPiNB();

    int* myStart = this->emissionParams->getStart();
    int myD = this->emissionParams->getD();
    int* myT = this->emissionParams->getT();
    int lower_n = 0;
    int upper_n = this->emissionParams->getNsample();

    if(currN != -1)
    {
        lower_n = currN;
        upper_n = currN+1;
    }
    int allT = 0;
    for (n = lower_n; n < upper_n; n++)
    {
        allT += myT[n];
    }
    int counter=0;
    int t;
    int obs_d = myStart[0];
    PROTECT(myd=NEW_INTEGER(1));
    INTEGER_POINTER(myd)[0] = obs_d+1;

    double** myGammaAux = this->emissionParams->getGammaAux();

/*PROTECT(myObs=NEW_NUMERIC(allT*2));
counter=0;
for (n = lower_n; n < upper_n; n++) {
    for(t=0; t<myT[n]*2; t++) {
        if(state2flag[statecouple] == 1) { // forward state
            if(t<myT[n]) {
                NUMERIC_POINTER(myObs)[counter] = observations[n][t][obs_d];
            }
            else {
                NUMERIC_POINTER(myObs)[counter] = observations[n][t-myT[n]][revop[obs_d]];
            }
}
else if(state2flag[statecouple] == -1) { // reverse state
if(t<myT[n]) {
NUMERIC_POINTER(myObs)[counter] = observations[n][t][revop[obs_d]];
}
else {
NUMERIC_POINTER(myObs)[counter] = observations[n][t-myT[n]][obs_d];
}
}
else { // undirected state
if(t<myT[n]) {
NUMERIC_POINTER(myObs)[counter] = observations[n][t][obs_d];
}
else {
NUMERIC_POINTER(myObs)[counter] = observations[n][t-myT[n]][obs_d];
}
}
counter++;
}
}*/
/*
    PROTECT(myGammas=NEW_NUMERIC(allT*2));
    counter = 0;
    for (n = lower_n; n < upper_n; n++) {
        for(t=0; t<myT[n]*2; t++) {
            if(state2flag[statecouple] == 1) { // forward state
                if(t<myT[n]) {
                    NUMERIC_POINTER(myGammas)[counter] = myGammaAux[n][t];
                }
                else {
                    NUMERIC_POINTER(myGammas)[counter] = revGammaAux[n][t-myT[n]];
}
}
else if(state2flag[statecouple] == -1) { // reverse state
if(t<myT[n]) {
NUMERIC_POINTER(myGammas)[counter] = revGammaAux[n][t];
}
else {
NUMERIC_POINTER(myGammas)[counter] = myGammaAux[n][t-myT[n]];
}
}
else { // undirected state
if(t<myT[n]) {
NUMERIC_POINTER(myGammas)[counter] = myGammaAux[n][t];
}
else {
NUMERIC_POINTER(myGammas)[counter] = myGammaAux[n][t-myT[n]];
}
}
counter++;
}
}
*/

    PROTECT(myGammas=NEW_NUMERIC(allT*2));
    counter = 0;
    for (n = lower_n; n < upper_n; n++)
    {
        for(t=0; t<myT[n]; t++)
        {
            if(state2flag[statecouple] == 1)      // forward state
            {
                if(t<myT[n])
                {
                    NUMERIC_POINTER(myGammas)[counter] = myGammaAux[n][t];
                }
                else
                {
                    NUMERIC_POINTER(myGammas)[counter] = revGammaAux[n][t-myT[n]];
                }
            }
            else if(state2flag[statecouple] == -1)// reverse state
            {
                if(t<myT[n])
                {
                    NUMERIC_POINTER(myGammas)[counter] = revGammaAux[n][t];
                }
                else
                {
                    NUMERIC_POINTER(myGammas)[counter] = myGammaAux[n][t-myT[n]];
                }
            }
            else                                  // undirected state
            {
                if(t<myT[n])
                {
                    NUMERIC_POINTER(myGammas)[counter] = myGammaAux[n][t];
                }
                else
                {
                    NUMERIC_POINTER(myGammas)[counter] = myGammaAux[n][t-myT[n]];
                }
            }
            counter++;
        }
    }
    for (n = lower_n; n < upper_n; n++)
    {
        for(t=myT[n]; t<2*myT[n]; t++)
        {
            if(state2flag[statecouple] == 1)      // forward state
            {
                if(t<myT[n])
                {
                    NUMERIC_POINTER(myGammas)[counter] = myGammaAux[n][t];
                }
                else
                {
                    NUMERIC_POINTER(myGammas)[counter] = revGammaAux[n][t-myT[n]];
                }
            }
            else if(state2flag[statecouple] == -1)// reverse state
            {
                if(t<myT[n])
                {
                    NUMERIC_POINTER(myGammas)[counter] = revGammaAux[n][t];
                }
                else
                {
                    NUMERIC_POINTER(myGammas)[counter] = myGammaAux[n][t-myT[n]];
                }
            }
            else                                  // undirected state
            {
                if(t<myT[n])
                {
                    NUMERIC_POINTER(myGammas)[counter] = myGammaAux[n][t];
                }
                else
                {
                    NUMERIC_POINTER(myGammas)[counter] = myGammaAux[n][t-myT[n]];
                }
            }
            counter++;
        }
    }

    SEXP MUNB, SIZENB, PINB, CURRN, wname;
    PROTECT(MUNB=NEW_NUMERIC(1));
    NUMERIC_POINTER(MUNB)[0] = this->emissionParams->getMuNB();
    PROTECT(SIZENB=NEW_NUMERIC(1));
    NUMERIC_POINTER(SIZENB)[0] = this->emissionParams->getSizeNB();
    PROTECT(PINB=NEW_NUMERIC(1));
    NUMERIC_POINTER(PINB)[0] = this->emissionParams->getPiNB();

//	printf("%f \n\n",this->emissionParams->getMuNB());
    PROTECT(CURRN=NEW_NUMERIC(1));
    NUMERIC_POINTER(CURRN)[0] = currN;

    SEXP sizeF;
    PROTECT(sizeF = NEW_NUMERIC(this->emissionParams->getNsample()));
    for (n = lower_n; n < upper_n; n++)
    {
        NUMERIC_POINTER(sizeF)[n] = this->getParameter()->getSizeFactorNB()[n];
    }

    SEXP NCORES;
    PROTECT(NCORES=NEW_NUMERIC(1));
    NUMERIC_POINTER(NCORES)[0] = ncores;

    SEXP CURRSTATE;
    PROTECT(CURRSTATE=NEW_NUMERIC(1));
    NUMERIC_POINTER(CURRSTATE)[0] = this->getParameter()->getCurrState()+1;

    SEXP myPars;
    PROTECT(myPars=NEW_LIST(10));
    SET_ELEMENT(myPars, 0, MUNB);
    SET_ELEMENT(myPars, 1, SIZENB);
    SET_ELEMENT(myPars, 2, PINB);
    SET_ELEMENT(myPars, 3, myGammas);
//SET_ELEMENT(myPars, 4, myObs);
    SET_ELEMENT(myPars, 4, myd);
    SET_ELEMENT(myPars, 5, CURRN);
    SET_ELEMENT(myPars, 6, getListElement(this->emissionParams->getUniqueCountSplit(), "countSplit"));
    SET_ELEMENT(myPars, 7, emissionPrior);
    SET_ELEMENT(myPars, 8, NCORES);
    SET_ELEMENT(myPars, 9, CURRSTATE);

    PROTECT(wname = NEW_CHARACTER(10));
    SET_STRING_ELT(wname, 0, mkChar("mu"));
    SET_STRING_ELT(wname, 1, mkChar("size"));
    SET_STRING_ELT(wname, 2, mkChar("pi"));
    SET_STRING_ELT(wname, 3, mkChar("gamma"));
//SET_STRING_ELT(wname, 4, mkChar("obs"));
    SET_STRING_ELT(wname, 4, mkChar("d"));
    SET_STRING_ELT(wname, 5, mkChar("currN"));
    SET_STRING_ELT(wname, 6, mkChar("uniqueCountSplit"));
    SET_STRING_ELT(wname, 7, mkChar("sizeFactor"));
    SET_STRING_ELT(wname, 8, mkChar("ncores"));
    SET_STRING_ELT(wname, 9, mkChar("currstate"));

    SET_NAMES(myPars, wname);
//SEXP call, res;
//#pragma omp critical
//{
    SEXP call = PROTECT( lang2(getListElement(this->emissionParams->getUniqueCountSplit(), "optimFct"), myPars ) ) ;
    SEXP res = PROTECT( eval( call, R_GlobalEnv ) ) ;
//}

    double myMuNew = REAL(res)[0];
    double mySizeNew = REAL(res)[1];
    double myPiNew = REAL(res)[2];

    this->emissionParams->setMuNB(myMuNew);
    this->emissionParams->setSizeNB(mySizeNew);
    this->emissionParams->setPiNB(myPiNew);

    UNPROTECT(14);
    int i,j;
    if(observations != NULL)
    {
        double** upobs = this->getParameter()->getUniqueObsProb();
        int** ulens = this->getParameter()->getUniqueLens();
        int nsample =  this->getParameter()->getN();
        int myD =  this->getParameter()->getD();
        double* myval = (double*)malloc(sizeof(double)*1);

        for(n=lower_n; n<upper_n; n++)
        {
            for(j=0; j<ulens[n][0]; j++)
            {
                if(upobs[n][j] != -1)
                {
                    myval[0] = (double)j;
                    upobs[n][j] = this->calcEmissionProbability(myval, -1, n);
                }
            }
        }

        free(myval);
    }

}


double NegativeBinomial::Prior(SEXP hyperparams)
{
    int i, j;
    double out = 1;
    return out;

}


void NegativeBinomial::update(double ***observations, double* Pk,
int** isNaN, SEXP emissionPrior, int currN, int ncores)
{

    int n;
    SEXP myObs, pars, myGammas, myd;

//Rprintf("Getting ready to optimize!\n");
    PROTECT(pars=NEW_NUMERIC(3));
    NUMERIC_POINTER(pars)[0] = this->emissionParams->getMuNB();
    NUMERIC_POINTER(pars)[1] = this->emissionParams->getSizeNB();
    NUMERIC_POINTER(pars)[2] = this->emissionParams->getPiNB();

    int* myStart = this->emissionParams->getStart();
    int myD = this->emissionParams->getD();
    int* myT = this->emissionParams->getT();
    int lower_n = 0;
    int upper_n = this->emissionParams->getNsample();

    if(currN != -1)
    {
        lower_n = currN;
        upper_n = currN+1;
    }
    int allT = 0;
    for (n = lower_n; n < upper_n; n++)
    {
        allT += myT[n];
    }
    int counter=0;
    int t;
    int obs_d = myStart[0];
    PROTECT(myd=NEW_INTEGER(1));
    INTEGER_POINTER(myd)[0] = obs_d+1;

    double** myGammaAux = this->emissionParams->getGammaAux();

/*PROTECT(myObs=NEW_NUMERIC(allT));
counter=0;
for (n = lower_n; n < upper_n; n++) {
    for(t=0; t<myT[n]; t++) {
        NUMERIC_POINTER(myObs)[counter] = observations[n][t][obs_d];
        counter++;
    }
}*/

    PROTECT(myGammas=NEW_NUMERIC(allT));
    counter = 0;
    for (n = lower_n; n < upper_n; n++)
    {
        for(t=0; t<myT[n]; t++)
        {
            NUMERIC_POINTER(myGammas)[counter] = myGammaAux[n][t];
            counter++;
        }
    }

    SEXP MUNB, SIZENB, PINB, CURRN, wname;
    PROTECT(MUNB=NEW_NUMERIC(1));
    NUMERIC_POINTER(MUNB)[0] = this->emissionParams->getMuNB();
    PROTECT(SIZENB=NEW_NUMERIC(1));
    NUMERIC_POINTER(SIZENB)[0] = this->emissionParams->getSizeNB();
    PROTECT(PINB=NEW_NUMERIC(1));
    NUMERIC_POINTER(PINB)[0] = this->emissionParams->getPiNB();

    PROTECT(CURRN=NEW_NUMERIC(1));
    NUMERIC_POINTER(CURRN)[0] = currN;

    SEXP sizeF;
    PROTECT(sizeF = NEW_NUMERIC(this->emissionParams->getNsample()));
    for (n = lower_n; n < upper_n; n++)
    {
        NUMERIC_POINTER(sizeF)[n] = this->getParameter()->getSizeFactorNB()[n];
    }

    SEXP NCORES;
    PROTECT(NCORES=NEW_NUMERIC(1));
    NUMERIC_POINTER(NCORES)[0] = ncores;

    SEXP CURRSTATE;
    PROTECT(CURRSTATE=NEW_NUMERIC(1));
    NUMERIC_POINTER(CURRSTATE)[0] = this->getParameter()->getCurrState()+1;

    SEXP myPars;
    PROTECT(myPars=NEW_LIST(10));
    SET_ELEMENT(myPars, 0, MUNB);
    SET_ELEMENT(myPars, 1, SIZENB);
    SET_ELEMENT(myPars, 2, PINB);
    SET_ELEMENT(myPars, 3, myGammas);
//SET_ELEMENT(myPars, 4, myObs);
    SET_ELEMENT(myPars, 4, myd);
    SET_ELEMENT(myPars, 5, CURRN);
    SET_ELEMENT(myPars, 6, getListElement(this->emissionParams->getUniqueCountSplit(), "countSplit"));
    SET_ELEMENT(myPars, 7, emissionPrior);
    SET_ELEMENT(myPars, 8, NCORES);
    SET_ELEMENT(myPars, 9, CURRSTATE);

    PROTECT(wname = NEW_CHARACTER(10));
    SET_STRING_ELT(wname, 0, mkChar("mu"));
    SET_STRING_ELT(wname, 1, mkChar("size"));
    SET_STRING_ELT(wname, 2, mkChar("pi"));
    SET_STRING_ELT(wname, 3, mkChar("gamma"));
//SET_STRING_ELT(wname, 4, mkChar("obs"));
    SET_STRING_ELT(wname, 4, mkChar("d"));
    SET_STRING_ELT(wname, 5, mkChar("currN"));
    SET_STRING_ELT(wname, 6, mkChar("uniqueCountSplit"));
    SET_STRING_ELT(wname, 7, mkChar("sizeFactor"));
    SET_STRING_ELT(wname, 8, mkChar("ncores"));
    SET_STRING_ELT(wname, 9, mkChar("currstate"));

    SET_NAMES(myPars, wname);
//SEXP call, res;
//#pragma omp critical
//{
    SEXP call = PROTECT( lang2(getListElement(this->emissionParams->getUniqueCountSplit(), "optimFct"), myPars ) ) ;
    SEXP res = PROTECT( eval( call, R_GlobalEnv ) ) ;
//}

    double myMuNew = REAL(res)[0];
    double mySizeNew = REAL(res)[1];
    double myPiNew = REAL(res)[2];

    this->emissionParams->setMuNB(myMuNew);
    this->emissionParams->setSizeNB(mySizeNew);
    this->emissionParams->setPiNB(myPiNew);

    UNPROTECT(14);
    int i,j,d;
    if(observations != NULL)
    {
        double** upobs = this->getParameter()->getUniqueObsProb();
        int** ulens = this->getParameter()->getUniqueLens();
        int nsample =  this->getParameter()->getN();
        int myD =  this->getParameter()->getD();
        double* myval = (double*)malloc(sizeof(double)*1);

        for(n=lower_n; n<upper_n; n++)
        {
            for(j=0; j<ulens[n][0]; j++)
            {
                if(upobs[n][j] != -1)
                {
                    myval[0] = (double)j;
                    upobs[n][j] = this->calcEmissionProbability(myval, -1, n);
                }
            }
        }

        free(myval);
    }
}


void NegativeBinomial::setParsToTwin(EmissionFunction* myTwinEmission, int currN, double*** observations)
{
    this->emissionParams->setMuNB(myTwinEmission->getParameter()->getMuNB());
    this->emissionParams->setSizeNB(myTwinEmission->getParameter()->getSizeNB());
    this->emissionParams->setPiNB(myTwinEmission->getParameter()->getPiNB());
    int i,j,d,n;
    int lower_n = 0;
    int upper_n = this->emissionParams->getNsample();
    if(currN != -1)
    {
        lower_n = currN;
        upper_n = currN+1;
    }
    if(observations != NULL)
    {
        double** upobs = this->getParameter()->getUniqueObsProb();
        int** ulens = this->getParameter()->getUniqueLens();
        int nsample =  this->getParameter()->getN();
        int myD =  this->getParameter()->getD();
        double* myval = (double*)malloc(sizeof(double)*1);

        for(n=lower_n; n<upper_n; n++)
        {
            for(j=0; j<ulens[n][0]; j++)
            {
                if(upobs[n][j] != -1)
                {
                    myval[0] = (double)j;
                    upobs[n][j] = this->calcEmissionProbability(myval, -1, n);
                }
            }
        }

        free(myval);
    }
//Rprintf("no\n");
}
