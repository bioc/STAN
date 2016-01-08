#include "PoissonLogNormal.h"
#include <cmath>
#define OBSERVATIONPOS(index, start) (start[index])

PoissonLogNormal::PoissonLogNormal(
ParamContainerEmissions *emissionParams)
{

    this->initial=1;
    this->emissionParams = emissionParams;

    int mem = sizeof(double) * this->emissionParams->getD() * 4
        + sizeof(double*) * this->emissionParams->getD() * 2;

    if (DEBUG_MEMORY)
    {
        printf("new->PoissonLogNormal; (%d bytes) \n", mem);
    }

}


PoissonLogNormal::~PoissonLogNormal()
{

    int d1;
    if (DEBUG_MEMORY)
    {
        printf("delete->PoissonLogNormal; (%d bytes) ", 0);
    }
    delete this->emissionParams;
}


double PoissonLogNormal::calcEmissionProbability(double *obs, int isna, int currN)
{

    int i, obs_i, D;
    double probability = 1;
    int* myStart = this->emissionParams->getStart();
    D = this->emissionParams->getD();
//printf("sF=%f\n",this->emissionParams->getSizeFactorPoiLog()[currN]);
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

            SEXP call, res, ALLPARS;

            PROTECT(ALLPARS=NEW_NUMERIC(3));
            NUMERIC_POINTER(ALLPARS)[0] = obs[obs_i];
            NUMERIC_POINTER(ALLPARS)[1] = this->emissionParams->getMuPoiLog()-log(this->emissionParams->getSizeFactorPoiLog()[currN]);;
            NUMERIC_POINTER(ALLPARS)[2] = this->emissionParams->getSigmaPoiLog();;

            call = PROTECT( lang2(install("call_dpoilog"), ALLPARS) ) ;
            res = PROTECT( eval( call, R_GlobalEnv ) ) ;
            double prob = REAL(res)[0];

            UNPROTECT(3);
// -log(this->emissionParams->getSizeFactorPoiLog()[currN])
//	probability = poilog((int)obs[obs_i], this->emissionParams->getMuPoiLog()-log(this->emissionParams->getSizeFactorPoiLog()[currN]), this->emissionParams->getSigmaPoiLog()*this->emissionParams->getSigmaPoiLog());
//	printf("%f == %f\n", prob, probability);

//probability = poilog((int)obs[obs_i], this->emissionParams->getMuPoiLog(), this->emissionParams->getSigmaPoiLog());
            probability = prob;
        }
    }
    else
    {
        obs_i = myStart[0];
        int myIndex = (int)obs[obs_i];
        probability = this->getParameter()->getUniqueObsProb()[currN][myIndex];
//	printf("%f\n", probability);
    }
//probability = poilog((int)obs[obs_i], this->emissionParams->getMuPoiLog()-log(this->emissionParams->getSizeFactorPoiLog()[currN]), this->emissionParams->getSigmaPoiLog());
    if(probability < 0)
    {
        error("Negative probabilitiy in PoissonLogNormal!");
    }
    if (probability < 1e-100)
    {
//	printf("here!\n");
        probability = 1e-100;
    }

    return probability;
}


void PoissonLogNormal::updateAuxiliaries(double*** observations,
double** gamma, double* Pk, int* T, int n, int i, int** isNaN)
{
// TODO !!!
/*int t;
for (t = 0; t < T[n]; t++) {
    this->emissionParams->setGammaAux(gamma[t][i], n, t);
}*/
}


void PoissonLogNormal::updateAuxiliariesCoupled(double*** observations,
double** gamma, double* Pk, int* T, int n, int i, int statecouple,
int** isNaN)
{

// TODO!!
/*for (t = 0; t < T[n]; t++) {
    this->emissionParams->setGammaAuxSEXP(gamma[t][i], n, t);
}*/
}


void PoissonLogNormal::updateAuxiliariesCoupledRevop(double*** observations,
double** gamma, double* Pk, int* T, int n, int i, int statecouple,
int* state2flag, int* revop, int** isNaN)
{
/* TODO!!!
    for (t = 0; t < T[n]; t++) {
        this->emissionParams->setGammaAuxSEXP(gamma[t][i], n, t);
    }*/
}


void PoissonLogNormal::updateCoupledRevop(double ***observations,
double* Pk, int statecouple, int* state2flag, int* revop,
double** revGammaAux, int** isNaN, SEXP emissionPrior, int currN, int ncores)
{
    int n;
    SEXP myObs, pars, myGammas, myd;
//Rprintf("Getting ready to optimize!\n");

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
    NUMERIC_POINTER(MUNB)[0] = this->emissionParams->getMuPoiLog();
    PROTECT(SIZENB=NEW_NUMERIC(1));
    NUMERIC_POINTER(SIZENB)[0] = this->emissionParams->getSigmaPoiLog();

    PROTECT(CURRN=NEW_NUMERIC(1));
    NUMERIC_POINTER(CURRN)[0] = currN;

    SEXP NCORES;
    PROTECT(NCORES=NEW_NUMERIC(1));
    NUMERIC_POINTER(NCORES)[0] = ncores;

    SEXP CURRSTATE;
    PROTECT(CURRSTATE=NEW_NUMERIC(1));
    NUMERIC_POINTER(CURRSTATE)[0] = this->getParameter()->getCurrState()+1;

    SEXP myPars;
    PROTECT(myPars=NEW_LIST(9));
    SET_ELEMENT(myPars, 0, MUNB);
    SET_ELEMENT(myPars, 1, SIZENB);
    SET_ELEMENT(myPars, 2, myGammas);
    SET_ELEMENT(myPars, 3, myd);
    SET_ELEMENT(myPars, 4, CURRN);
    SET_ELEMENT(myPars, 5, getListElement(this->emissionParams->getUniqueCountSplit(), "countSplit"));
    SET_ELEMENT(myPars, 6, NCORES);
    SET_ELEMENT(myPars, 7, CURRSTATE);
    SET_ELEMENT(myPars, 8, emissionPrior);

    PROTECT(wname = NEW_CHARACTER(9));
    SET_STRING_ELT(wname, 0, mkChar("mu"));
    SET_STRING_ELT(wname, 1, mkChar("sigma"));
    SET_STRING_ELT(wname, 2, mkChar("gamma"));
    SET_STRING_ELT(wname, 3, mkChar("d"));
    SET_STRING_ELT(wname, 4, mkChar("currN"));
    SET_STRING_ELT(wname, 5, mkChar("uniqueCountSplit"));
    SET_STRING_ELT(wname, 6, mkChar("ncores"));
    SET_STRING_ELT(wname, 7, mkChar("currstate"));
    SET_STRING_ELT(wname, 8, mkChar("sizeFactor"));
    SET_NAMES(myPars, wname);

    SEXP call, res;
//	#pragma omp critical
//	{
    call = PROTECT( lang2(getListElement(this->emissionParams->getUniqueCountSplit(), "optimFct"), myPars ) ) ;
    res = PROTECT( eval( call, R_GlobalEnv ) ) ;
//	}

    double myMuNew = REAL(res)[0];
    double mySigmaNew = REAL(res)[1];
    this->emissionParams->setMuPoiLog(myMuNew);
    this->emissionParams->setSigmaPoiLog(mySigmaNew);

    UNPROTECT(11);
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


double PoissonLogNormal::Prior(SEXP hyperparams)
{
    int i, j;
    double out = 1;
    return out;

}


void PoissonLogNormal::update(double ***observations, double* Pk,
int** isNaN, SEXP emissionPrior, int currN, int ncores)
{

    int n;
    SEXP myObs, pars, myGammas, myd;

//Rprintf("Getting ready to optimize!\n");

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
    NUMERIC_POINTER(MUNB)[0] = this->emissionParams->getMuPoiLog();
    PROTECT(SIZENB=NEW_NUMERIC(1));
    NUMERIC_POINTER(SIZENB)[0] = this->emissionParams->getSigmaPoiLog();

    PROTECT(CURRN=NEW_NUMERIC(1));
    NUMERIC_POINTER(CURRN)[0] = currN;

    SEXP NCORES;
    PROTECT(NCORES=NEW_NUMERIC(1));
    NUMERIC_POINTER(NCORES)[0] = ncores;

    SEXP CURRSTATE;
    PROTECT(CURRSTATE=NEW_NUMERIC(1));
    NUMERIC_POINTER(CURRSTATE)[0] = this->getParameter()->getCurrState()+1;

    SEXP myPars;
    PROTECT(myPars=NEW_LIST(9));
    SET_ELEMENT(myPars, 0, MUNB);
    SET_ELEMENT(myPars, 1, SIZENB);
    SET_ELEMENT(myPars, 2, myGammas);
    SET_ELEMENT(myPars, 3, myd);
    SET_ELEMENT(myPars, 4, CURRN);
    SET_ELEMENT(myPars, 5, getListElement(this->emissionParams->getUniqueCountSplit(), "countSplit"));
    SET_ELEMENT(myPars, 6, NCORES);
    SET_ELEMENT(myPars, 7, CURRSTATE);
    SET_ELEMENT(myPars, 8, emissionPrior);

    PROTECT(wname = NEW_CHARACTER(9));
    SET_STRING_ELT(wname, 0, mkChar("mu"));
    SET_STRING_ELT(wname, 1, mkChar("sigma"));
    SET_STRING_ELT(wname, 2, mkChar("gamma"));
    SET_STRING_ELT(wname, 3, mkChar("d"));
    SET_STRING_ELT(wname, 4, mkChar("currN"));
    SET_STRING_ELT(wname, 5, mkChar("uniqueCountSplit"));
    SET_STRING_ELT(wname, 6, mkChar("ncores"));
    SET_STRING_ELT(wname, 7, mkChar("currstate"));
    SET_STRING_ELT(wname, 8, mkChar("sizeFactor"));
    SET_NAMES(myPars, wname);

    SEXP call, res;
//	#pragma omp critical
//	{
    call = PROTECT( lang2(getListElement(this->emissionParams->getUniqueCountSplit(), "optimFct"), myPars ) ) ;
    res = PROTECT( eval( call, R_GlobalEnv ) ) ;
//	}

    double myMuNew = REAL(res)[0];
    double mySigmaNew = REAL(res)[1];
    this->emissionParams->setMuPoiLog(myMuNew);
    this->emissionParams->setSigmaPoiLog(mySigmaNew);

    UNPROTECT(11);
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


void PoissonLogNormal::setParsToTwin(EmissionFunction* myTwinEmission, int currN, double*** observations)
{
    this->emissionParams->setMuPoiLog(myTwinEmission->getParameter()->getMuPoiLog());
    this->emissionParams->setSigmaPoiLog(myTwinEmission->getParameter()->getSigmaPoiLog());

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
}
