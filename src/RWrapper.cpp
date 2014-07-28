/**
 * @file    RWrapper.cpp
 * @author  Benedikt Zacher, AG Tresch, Gene Center Munich (zacher@lmb.uni-muenchen.de)
 * @version 0.99.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 *
 * This file implements the C wrapper function for interfacing C++ code from R.
 */

#include "matUtils.h"
#include "MemoryAllocation.h"
#include "ParamContainerEmissions.h"
#include "EmissionFactory.h"
#include "TransitionMatrix.h"
#include "InitialProbability.h"
#include "HMM.h"
#include "DebugConstants.h"
#include <Rembedded.h>
#include <Rdefines.h>

extern "C"
{
#include "RWrapper.h"
#include "MemoryAllocation.h"
#include "ParamContainerEmissions.h"
#include <Rembedded.h>
#include <Rdefines.h>
#include <list>
#define ARRAYSIZE(a) (sizeof(a) / sizeof(a[0]))
#define COLSIZE(a) ARRAYSIZE(a[0])

    R_CallMethodDef callMethods[] =
    {
        {"RHMMVITERBI", (DL_FUNC)&RHMMVITERBI, 9},
        {"RHMMFit", (DL_FUNC)&RHMMFit, 22},
        {"RGETPOSTERIOR", (DL_FUNC)&RGETPOSTERIOR, 12},
        {NULL, NULL, 0}
    };

    void R_init_STAN(DllInfo *info)
    {
        R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    }

    EmissionFunction** RGETMULTGAUSS(SEXP sexpmu, SEXP sexpsigma, int D, SEXP sexpk,int* start)
    {
        int parallel = 0;

        int i,j,k;

        int K = INTEGER(sexpk)[0];

        EmissionFactory* factory = createEmissionFactory(MULTIVARIATEGAUSSIAN);
        EmissionFunction **HMMEmissionFunctions = allocateEmissionFunctionVector(K);

        for(k=0; k<K; k++)
        {
            if(DEBUG_MEMORY)
            {
                Rprintf("\tstate=%d: ", k);
            }
            double **mu = allocateNumericMatrix(D, 1);
            for(i=0; i<D; i++)
            {
                mu[i][0] = REAL(VECTOR_ELT(sexpmu, k))[i];
            }

            double **sigma = allocateNumericMatrix(D, D);
            for(i=0; i<D; i++)
            {
                for(j=0; j<D; j++)
                {
                    sigma[i][j] = REAL(coerceVector(VECTOR_ELT(sexpsigma, k), REALSXP))[i+D*j];
                }
            }

            HMMEmissionFunctions[k] = factory->createEmissionFunction(new ParamContainerEmissions(mu, sigma, 0, D,start), parallel);
        }

        delete factory;

        return HMMEmissionFunctions;
    }

    EmissionFunction** RGETBERNOULLI2(SEXP sexpbernoullip, int D, SEXP sexpk, int* start)
    {
        int i,j,k;
        int K = INTEGER(sexpk)[0];

        EmissionFactory* factory = createEmissionFactory(BERNOULLI);
        EmissionFunction **HMMEmissionFunctions = allocateEmissionFunctionVector(K);
        for(k=0; k<K; k++)
        {
            double p1 = REAL(coerceVector(VECTOR_ELT(sexpbernoullip, k), REALSXP))[0];
            // create new constructor in ParamContainerEmissions (p* and dim(int) )
            HMMEmissionFunctions[k] = factory->createEmissionFunction(new ParamContainerEmissions(p1, D, start), 0);
        }

        delete factory;

        return HMMEmissionFunctions;
    }

    EmissionFunction** RGETEMISSION(SEXP sexpparameters, int D, SEXP sexpk, int* start, const char* type)
    {
        EmissionFunction **HMMEmissionFunctions;
        const char* gaussian = "Gaussian";
        const char* bernoulli = "Bernoulli";

        SEXP t=R_NilValue;
        PROTECT(t = GET_SLOT(sexpparameters, install("parameters")));
        GET_SLOT(sexpparameters, install("parameters"));
        // TODO: heir muss ich den parameter slot rausziehen
        if(strcmp(type,bernoulli) == 0)
        {
            HMMEmissionFunctions = RGETBERNOULLI2(getListElement(t, "p"), D, sexpk, start);
        }
        else if(strcmp(type,gaussian) == 0)
        {
            HMMEmissionFunctions = RGETMULTGAUSS(getListElement(t, "mean"), getListElement(t, "cov"), D, sexpk, start);
        }
        UNPROTECT(1);

        return HMMEmissionFunctions;
    }


    EmissionFunction** RGETBERNOULLI(SEXP sexpbernoullip, int D, SEXP sexpk, int* start, int row)
    {
        int i,j,k;
        int K = INTEGER(sexpk)[0];

        EmissionFactory* factory = createEmissionFactory(BERNOULLI);
        EmissionFunction **HMMEmissionFunctions = allocateEmissionFunctionVector(K);
        SEXP innerP, inP;

        for(k=0; k<K; k++)
        {
            //for p's
            double p1;
            if (row == -1)
            {
                innerP = getListElement(sexpbernoullip, "p");
                p1 = REAL(coerceVector(VECTOR_ELT(innerP, k), REALSXP))[0];
            }
            else
            {
                inP = VECTOR_ELT(sexpbernoullip, row);
                p1= REAL(coerceVector(VECTOR_ELT(inP, k), REALSXP))[0];
            }

            // create new constructor in ParamContainerEmissions (p* and dim(int) ) //0  war parallel param
            HMMEmissionFunctions[k] = factory->createEmissionFunction(new ParamContainerEmissions(p1, D, start), 0);
        }

        delete factory;

        return HMMEmissionFunctions;
    }

    SEXP RPREPAREGAUSSPAR(EmissionFunction** myEmissions, int K, int useNames)
    {
        int k,i,j;
        int D =  myEmissions[0]->getParameter()->getD();

        SEXP muFit, sigmaFit, inverseSigmaFit, emissionParam, wnames, currMU, currSIGMA, currINVSIGMA;
        PROTECT(emissionParam = NEW_LIST(3));

        PROTECT(muFit = NEW_LIST(K));
        for(k=0; k<K; k++)
        {
            PROTECT(currMU = NEW_NUMERIC(D));
            for(i=0; i<D; i++)
            {
                NUMERIC_POINTER(currMU)[i] = myEmissions[k]->getParameter()->getGaussianMU()[i][0];

            }
            SET_ELEMENT(muFit, k, currMU);
        }

        PROTECT(sigmaFit = NEW_LIST(K));
        for(k=0; k<K; k++)
        {
            PROTECT(currSIGMA = NEW_NUMERIC(D*D));
            for(i=0; i<D; i++)
            {
                for(j=0; j<D; j++)
                {
                    NUMERIC_POINTER(currSIGMA)[j+D*i] = myEmissions[k]->getParameter()->getGaussianSIGMA()[i][j];
                }
            }
            SET_ELEMENT(sigmaFit, k, currSIGMA);
        }

        PROTECT(inverseSigmaFit = NEW_LIST(K));
        for(k=0; k<K; k++)
        {
            PROTECT(currINVSIGMA = NEW_NUMERIC(D*D));
            for(i=0; i<D; i++)
            {
                for(j=0; j<D; j++)
                {
                    NUMERIC_POINTER(currINVSIGMA)[j+D*i] = myEmissions[k]->getParameter()->getGaussianINVSIGMA()[i][j];
                }
            }
            SET_ELEMENT(inverseSigmaFit, k, currINVSIGMA);
        }

        if(useNames)
        {
            PROTECT(wnames = NEW_CHARACTER(3));
            SET_STRING_ELT(wnames, 0, mkChar("mean"));
            SET_STRING_ELT(wnames, 1, mkChar("cov"));
            SET_STRING_ELT(wnames, 2, mkChar("invsigma"));
            SET_NAMES(emissionParam, wnames);
            UNPROTECT(1);
        }

        SET_ELEMENT(emissionParam, 0, muFit);
        SET_ELEMENT(emissionParam, 1, sigmaFit);
        SET_ELEMENT(emissionParam, 2, inverseSigmaFit);
        UNPROTECT(4+3*K);

        return emissionParam;
    }

    SEXP RPREPAREBERNOULLIPAR(EmissionFunction** myEmissions, int K)
    {
        int k,i,j;
        int D = myEmissions[0]->getParameter()->getD();

        SEXP emissionParam, currP, pFit, wname;
        PROTECT(emissionParam = NEW_LIST(1));

        PROTECT(pFit = NEW_LIST(K));
        for(k=0; k<K; k++)
        {
            PROTECT(currP = NEW_NUMERIC(D));
            for(i=0; i<D; i++)
            {
                NUMERIC_POINTER(currP)[i] = myEmissions[k]->getParameter()->getBernoulliP();
            }
            SET_ELEMENT(pFit, k, currP);
        }

        PROTECT(wname = NEW_CHARACTER(1));
        SET_STRING_ELT(wname, 0, mkChar("p"));
        SET_NAMES(emissionParam, wname);
        UNPROTECT(1);

        SET_ELEMENT(emissionParam, 0, pFit);
        UNPROTECT(2+1*K);

        return emissionParam;

    }

    SEXP RPREPAREBERNOULLIPAR2(EmissionFunction** myEmissions, int K, int useNames)
    {
        int k,d,j;
        int D = myEmissions[0]->getParameter()->getD();

        SEXP p, currState, wname;

        PROTECT(p = NEW_LIST(K));
        for(k=0; k<K; k++)
        {
            PROTECT(currState = NEW_NUMERIC(D));
            for(d=0; d<D; d++)
            {
                NUMERIC_POINTER(currState)[d] = myEmissions[k]->getParameter()->getBernoulliP();
            }
            SET_ELEMENT(p, k, currState);
        }

        if(useNames)
        {
            PROTECT(wname = NEW_CHARACTER(1));
            SET_STRING_ELT(wname, 0, mkChar("p"));
            SET_NAMES(p, wname);
            UNPROTECT(1);
        }

        UNPROTECT(K+1);

        return p;

    }

    SEXP RPREPAREEMISSIONPAR(EmissionFunction** myEmissions, int K, const char* type, int useNames)
    {
        if(strcmp(type,"Gaussian") == 0)
        {
            return RPREPAREGAUSSPAR(myEmissions, K, useNames);
        }
        else if(strcmp(type,"Bernoulli") == 0)
        {
            return RPREPAREBERNOULLIPAR2(myEmissions, K, useNames);
        }
    }

    SEXP RPREPAREJOINTLYINDEPENDENTPAR(EmissionFunction** myEmissions, int K, SEXP types)
    {
        SEXP emissionList, stateList;
        int d,k;

        PROTECT(emissionList = NEW_LIST(LENGTH(types)));
        for(d=0; d<LENGTH(types); d++)
        {
            PROTECT(stateList = NEW_LIST(K));
            SET_ELEMENT(emissionList, d, stateList);
        }

        std::list<EmissionFunction*> listEF;
        for (k = 0; k < K; k++)
        {
            std::list<EmissionFunction*>::iterator pos;
            listEF = myEmissions[k]->getEmissionFunctionList();
            d = 0;
            for (pos = listEF.begin(); pos!= listEF.end(); pos++)
            {
                const char* type = CHAR(STRING_ELT(types, d));
                EmissionFunction** currEmission = (EmissionFunction**)malloc(sizeof(EmissionFunction*));
                currEmission[0] = (*pos);
                SET_ELEMENT(VECTOR_ELT(emissionList, d), k, RPREPAREEMISSIONPAR(currEmission, 1, type, 0));
                d++;
            }
        }

        SEXP output;
        PROTECT(output = NEW_LIST(2));
        SET_ELEMENT(output, 0, emissionList);
        SET_ELEMENT(output, 1, types);

        SEXP wname;
        PROTECT(wname = NEW_CHARACTER(2));
        SET_STRING_ELT(wname, 0, mkChar("emissions"));
        SET_STRING_ELT(wname, 1, mkChar("types"));
        SET_NAMES(output, wname);

        UNPROTECT(3+LENGTH(types));

        return output;
    }

    /*
     * JOINTLY INDEPENDENT
     *
     */
    EmissionFunction** createJointlyIndependent(std::list<EmissionFunction**> bernoulligauss, int D, SEXP sexpk, int* T, int nsample)
    {
        int i,j,k;
        int K = INTEGER(sexpk)[0];

        EmissionFunction **HMMEmissionFunctions = allocateEmissionFunctionVector(K);
        EmissionFactory* factory = createEmissionFactory(JOINTLYINDEPENDENT);

        for(k=0; k<K; k++)
        {
            if(DEBUG_MEMORY)
            {
                Rprintf("\tstate=%d: ", k);
            }
            std::list<EmissionFunction**>::iterator pos;
            std::list<EmissionFunction*> combinedPerState;
            for (pos = bernoulligauss.begin(); pos!=bernoulligauss.end(); pos++)
            {
                EmissionFunction* f = (*pos)[k];
                combinedPerState.push_back(f);
            }

            HMMEmissionFunctions[k] = factory->createEmissionFunctionMixed(combinedPerState, new ParamContainerEmissions(D));
            HMMEmissionFunctions[k]->getParameter()->setDataVars(nsample, T);

            list<EmissionFunction*>::iterator subEmission;
            std::list<EmissionFunction*> listEF = HMMEmissionFunctions[k]->getEmissionFunctionList();
            // set gamma aux of "sub-emissions" to gamma aux of wrapper emission
            for (subEmission = listEF.begin(); subEmission != listEF.end(); subEmission++)
            {
                double **myGammaAux = HMMEmissionFunctions[k]->getParameter()->getGammaAux();
                ((*subEmission)->getParameter())->setDataVars(myGammaAux);
            }
        }

        delete factory;
        return HMMEmissionFunctions;
    }

    SEXP RPREPAREMIXEDBERNOULLIGAUSSPAR(EmissionFunction** myEmissions, int K, SEXP types, SEXP order)
    {
        int k,i,j,z, bernoulli, gaussian, unprotectVal, gaussianLength, bernoulliLength, stackCounter;
        bernoulli = 0;
        gaussian = 0;
        gaussianLength = 0;
        z= 0;
        unprotectVal = 0;
        bernoulliLength= 0;
        stackCounter=0;
        SEXP emissionParam, currP, pFit, wname, pvalues, pname, Fitorder,Ordernames;
        SEXP muFit, sigmaFit, Fittypes, inverseSigmaFit, wnames, currMU,currOrder, currSIGMA, currINVSIGMA, pFitAll, pnames;

        for (j = 0; j< length(types); j++)
        {
            const char* typ = CHAR(STRING_ELT(types, j));
            if (strcmp(typ,"Gaussian") == 0 && gaussian == 0)
            {
                PROTECT(muFit = NEW_LIST(K));
                PROTECT(sigmaFit = NEW_LIST(K));
                PROTECT(inverseSigmaFit = NEW_LIST(K));
                gaussian = 1;
                stackCounter=stackCounter+3;

            }
            if (strcmp(typ,"Gaussian") == 0 )
            {
                gaussianLength++;
            }
            if (strcmp(typ,"Bernoulli") == 0 && bernoulli == 0)
            {
                bernoulli = 1;
            }
            if (strcmp(typ,"Bernoulli") == 0)
            {
                bernoulliLength++;
            }

        }
        if (bernoulliLength > 0 )
        {
            PROTECT(pFitAll=NEW_LIST(bernoulliLength));
            PROTECT(pnames = NEW_CHARACTER(bernoulliLength));
            stackCounter=stackCounter+2;
        }
        else
        {
            bernoulliLength = 1;
        }

        if (bernoulli ==1 && gaussian==1)
        {
            PROTECT(emissionParam = NEW_LIST(6));
            stackCounter=stackCounter+1;
        }
        if (bernoulli ==1 && gaussian == 0)
        {
            PROTECT(emissionParam = NEW_LIST(3));
            stackCounter=stackCounter+1;
        }
        if (bernoulli ==0 && gaussian == 1)
        {
            PROTECT(emissionParam = NEW_LIST(5));
            stackCounter=stackCounter+1;
        }

        std::list<EmissionFunction*> listEF;
        int sizeLEF, b, l;
        l= 0;

        for (b = 0; b< bernoulliLength; b++)
        {

            if (bernoulli == 1 )
            {
                SET_STRING_ELT(pnames, b, mkChar("p"));
                PROTECT(pFit= NEW_LIST(K));
                stackCounter=stackCounter+1;

            }

            for (k = 0; k < K; k++)
            {
                //extract EmissionFunctions for state k
                std::list<EmissionFunction*>::iterator pos;
                listEF =myEmissions[k]->getEmissionFunctionList();
                int t = 0;
                for (pos = listEF.begin(); pos!= listEF.end(); pos++)
                {
                    const char* typ = CHAR(STRING_ELT(types, t));
                    t++;
                    int D = (*pos)->getParameter()->getD();

                    if (strcmp(typ,"Gaussian") == 0)
                    {

                        PROTECT(currMU = NEW_NUMERIC((*pos)->getParameter()->getD()));
                        for (i = 0; i < (*pos)->getParameter()->getD(); i++)
                        {
                            NUMERIC_POINTER(currMU)[i] =(*pos)->getParameter()->getGaussianMU()[i][0];

                        }
                        SET_ELEMENT(muFit, k, currMU);
                        PROTECT(currSIGMA = NEW_NUMERIC(D * D));
                        for (i = 0; i < D; i++)
                        {
                            for (j = 0; j < D; j++)
                            {
                                NUMERIC_POINTER(currSIGMA)[j + D * i] =
                                    (*pos)->getParameter()->getGaussianSIGMA()[i][j];
                            }
                        }
                        SET_ELEMENT(sigmaFit, k, currSIGMA);

                        PROTECT(currINVSIGMA = NEW_NUMERIC(D * D));
                        for (i = 0; i < D; i++)
                        {
                            for (j = 0; j < D; j++)
                            {
                                NUMERIC_POINTER(currINVSIGMA)[j + D * i] =
                                    (*pos)->getParameter()->getGaussianINVSIGMA()[i][j];
                            }
                        }
                        SET_ELEMENT(inverseSigmaFit, k, currINVSIGMA);
                        stackCounter=stackCounter+3;

                    }
                }
                sizeLEF = t;
                t = 0;
                int z;
                z = 0;

                for (pos = listEF.begin(); pos!= listEF.end(); pos++)
                {
                    const char* typ = CHAR(STRING_ELT(types, t));
                    t++;
                    if (strcmp(typ,"Bernoulli") == 0 )
                    {
                        if ( z == b )
                        {
                            PROTECT(currP = NEW_NUMERIC(1));
                            NUMERIC_POINTER(currP)[0] =(*pos)->getParameter()->getBernoulliP();
                            SET_ELEMENT(pFit, k, currP);
                            stackCounter=stackCounter+1;
                        }
                        z++;
                    }                             
                    // end of one Bernoulli, add to pFitAll
                }
            }
            // end of k States

            l++;
            if (bernoulli == 1 )
            {
                SET_ELEMENT(pFitAll, b, pFit);
            }

        }                                         
        // end of all Bernoullis

        if (bernoulli == 1 )
        {
            SET_NAMES(pFitAll, pnames);
        }

        z=0;
        if (bernoulli == 1 && gaussian==1)
        {
            sizeLEF = bernoulliLength + gaussianLength;
        }

        if (bernoulli== 0 && gaussian==1 )
        {
            sizeLEF = gaussianLength;
        }
        if (bernoulli==1 && gaussian==0 )
        {
            sizeLEF = bernoulliLength;
        }

        int a, h;
        a =0;
        PROTECT(Fitorder=NEW_LIST(sizeLEF));
        PROTECT(Ordernames=NEW_CHARACTER(sizeLEF));
        stackCounter=stackCounter+2;
        for (a= 0; a< sizeLEF; a++ )
        {
            int ordersize;
            SEXP vect = VECTOR_ELT(order, a);
            ordersize = length(vect);
            PROTECT(currOrder = NEW_NUMERIC(ordersize));
            const char* typ = CHAR(STRING_ELT(types, a));
            for (h = 0; h < ordersize; h++)
            {
                NUMERIC_POINTER(currOrder)[h] = INTEGER(vect)[h];

            }
            SET_STRING_ELT(Ordernames, a, mkChar(typ));
            SET_ELEMENT(Fitorder, a, currOrder);
            UNPROTECT(1);

        }
        SET_NAMES(Fitorder, Ordernames);

        // states are over
        if (bernoulli ==1 && gaussian==1)
        {
            PROTECT(wnames = NEW_CHARACTER(6));
            SET_STRING_ELT(wnames, 0, mkChar("mean"));
            SET_STRING_ELT(wnames, 1, mkChar("cov"));
            SET_STRING_ELT(wnames, 2, mkChar("invsigma"));
            SET_STRING_ELT(wnames, 3, mkChar("Bernoulli"));
            SET_STRING_ELT(wnames, 4, mkChar("order"));
            SET_STRING_ELT(wnames, 5, mkChar("types"));
            SET_ELEMENT(emissionParam, 0, muFit);
            SET_ELEMENT(emissionParam, 1, sigmaFit);
            SET_ELEMENT(emissionParam, 2, inverseSigmaFit);
            SET_ELEMENT(emissionParam, 3, pFitAll);
            SET_ELEMENT(emissionParam, 4, Fitorder);
            SET_ELEMENT(emissionParam, 5, Ordernames);
        }
        if (bernoulli ==1 && gaussian == 0)
        {
            PROTECT(wnames = NEW_CHARACTER(3));
            SET_STRING_ELT(wnames, 0, mkChar("Bernoulli"));
            SET_STRING_ELT(wnames, 1, mkChar("order"));
            SET_STRING_ELT(wnames, 2, mkChar("types"));
            SET_ELEMENT(emissionParam, 0, pFitAll);
            SET_ELEMENT(emissionParam, 1, Fitorder);
            SET_ELEMENT(emissionParam, 2, Ordernames);
        }
        if (bernoulli ==0 && gaussian == 1)
        {
            PROTECT(wnames = NEW_CHARACTER(5));
            SET_STRING_ELT(wnames, 0, mkChar("mean"));
            SET_STRING_ELT(wnames, 1, mkChar("cov"));
            SET_STRING_ELT(wnames, 2, mkChar("invsigma"));
            SET_STRING_ELT(wnames, 3, mkChar("order"));
            SET_STRING_ELT(wnames, 4, mkChar("types"));
            SET_ELEMENT(emissionParam, 0, muFit);
            SET_ELEMENT(emissionParam, 1, sigmaFit);
            SET_ELEMENT(emissionParam, 2, inverseSigmaFit);
            SET_ELEMENT(emissionParam, 3, Fitorder);
            SET_ELEMENT(emissionParam, 4, Ordernames);
        }

        SET_NAMES(emissionParam, wnames);

        UNPROTECT(stackCounter+1);

        return emissionParam;

    }

    double*** RGETOBS(SEXP sexpobs, int* T, int nsample, int D)
    {
        int t,d,n;
        double ***observations = NULL;

        if(nsample > 0)
        {
            int mem = 0;
            observations = (double***)malloc(sizeof(double**)*nsample);
            mem += sizeof(double**)*nsample;
            for(n=0; n<nsample; n++)
            {
                observations[n] = (double**)malloc(sizeof(double*)*T[n]);
                mem += sizeof(double*)*T[n];
                for(t=0; t<T[n]; t++)
                {
                    observations[n][t] = (double*)malloc(sizeof(double)*D);
                    mem += sizeof(double)*D;
                    for(d=0; d<D; d++)
                    {
                        observations[n][t][d] = REAL(coerceVector(VECTOR_ELT(sexpobs, n), REALSXP))[t+T[n]*d];
                    }
                }
            }
            if(DEBUG_MEMORY)
            {
                Rprintf("new->observation matrix; (%d bytes)\n", mem);
            }
        }
        return observations;
    }

    TransitionMatrix* RGETTRANSMAT(SEXP sexpA, int K)
    {
        int i,j;

        sexpA = coerceVector(sexpA, REALSXP);
        double **transProb = (double**)malloc(sizeof(double*)*K);

        for(i=0; i<K; i++)
        {
            transProb[i] = (double*)malloc(sizeof(double)*K);
            for(j=0; j<K; j++)
            {
                transProb[i][j] = REAL(sexpA)[i+K*j];
            }
        }

        return new TransitionMatrix(transProb, K);
    }

    InitialProbability* RGETINITPROB(SEXP sexppi, int K)
    {
        int i;
        double *pi = (double*)malloc(sizeof(double)*K);
        for(i=0; i<K; i++)
        {
            pi[i] = REAL(sexppi)[i];
        }
        return new InitialProbability(pi, K);
    }

    SEXP RPREPAREPI(InitialProbability* myInitProb)
    {
        SEXP piFit;
        int i;

        PROTECT(piFit = NEW_NUMERIC(myInitProb->getK()));
        for(i=0; i<myInitProb->getK(); i++)
        {
            NUMERIC_POINTER(piFit)[i] = myInitProb->getInitialProb()[i];
        }
        UNPROTECT(1);
        return piFit;
    }

    SEXP RPREPARETRANSMAT(TransitionMatrix* transMat)
    {
        SEXP AFit;
        int i,j;
        int K = transMat->getK();
        PROTECT(AFit = NEW_NUMERIC(K*K));
        for(i=0; i<K; i++)
        {
            for(j=0; j<K; j++)
            {
                NUMERIC_POINTER(AFit)[j+K*i] = transMat->getTransMat()[i][j];
            }
        }
        UNPROTECT(1);
        return AFit;
    }

    HMM* createHMM(int parallel, int K, InitialProbability* initProb, TransitionMatrix* transMat, EmissionFunction** myEmissions)
    {
        if(parallel == 0)
        {
            return new HMM(K, initProb, transMat, myEmissions);
        }
    }

    void RGETFLAGS(SEXP sexpflags, SEXP sexpstate2flag, int*** flags, int** state2flag, int nsample, int* T, int K)
    {
        int n,t,i;

        if(LENGTH(sexpflags) != 0)
        {
            *flags = (int**)malloc(sizeof(int*)*nsample);
            for(n=0; n<nsample; n++)
            {
                (*flags)[n] = (int*)malloc(sizeof(int)*T[n]);
                for(t=0; t<T[n]; t++)
                {
                    (*flags)[n][t] = INTEGER(VECTOR_ELT(sexpflags, n))[t];
                }
            }
        }

        if(LENGTH(sexpstate2flag) != 0)
        {
            *state2flag = (int*)malloc(sizeof(int)*K);
            for(i=0; i<K; i++)
            {
                (*state2flag)[i] = INTEGER(sexpstate2flag)[i];
            }
        }
    }

    void RFREEFLAGS(SEXP sexpflags, SEXP sexpstate2flag, int** flags, int* state2flag, int nsample)
    {
        int n,t,i;

        if(LENGTH(sexpflags) != 0)
        {
            for(n=0; n<nsample; n++)
            {
                free(flags);
            }
            free(flags);

        }

        if(LENGTH(sexpstate2flag) != 0)
        {
            free(state2flag);
        }

    }

    void RGETCOUPLES(SEXP sexpcouples, int** couples, int K)
    {
        int i;

        if(LENGTH(sexpcouples) != 0)
        {
            *couples = (int*)malloc(sizeof(int)*K);
            for(i=0; i<K; i++)
            {
                (*couples)[i] = INTEGER(sexpcouples)[i];
            }
        }
    }

    int** whichNaN(double*** obs, int nsample, int* T, int D)
    {

        int n,t,d;
        int** myNaN = NULL;

        if(nsample > 0)
        {
            myNaN = (int**)malloc(sizeof(int*)*nsample);

            for(n=0; n<nsample; n++)
            {
                myNaN[n] = (int*)malloc(sizeof(int)*T[n]);
                for(t=0; t<T[n]; t++)
                {
                    myNaN[n][t] = 0;
                    for(d=0; d<D; d++)
                    {
                        if(obs[n][t][d] != obs[n][t][d])
                        {
                            myNaN[n][t] = 1;
                        }
                    }
                }
            }
        }
        return myNaN;
    }

    EmissionFunction** getEmission(const char* type, SEXP sexpemission, SEXP sexpk, int* start_d, int nsample, int* T, int K, int D)
    {
        int i;
        EmissionFunction** myEmissions = NULL;
        const char* gauss = "Gaussian";
        const char* independent = "Independent";
        const char* jointlyindependent = "JointlyIndependent";

        if(strcmp(type, gauss) == 0)
        {
            //parallel
            myEmissions = RGETMULTGAUSS(getListElement(sexpemission, "mean"), getListElement(sexpemission, "cov"), D, sexpk, start_d);
            for(i=0; i<K; i++)
            {
                myEmissions[i]->getParameter()->setDataVars(nsample, T);
            }
        }
        else if(strcmp(type, jointlyindependent) == 0)
        {
            int nemissions = LENGTH(getListElement(sexpemission, "emissions"));
            list<EmissionFunction**> combinedDistributions;

            SEXP myInpedendentEmissions = getListElement(sexpemission, "emissions");
            int currEmission;
            for (currEmission = 0; currEmission < nemissions; currEmission++)
            {
                int* Ds = NULL;
                int ordersize;
                SEXP currDims = getListElement(sexpemission, "emissionDim");
                SEXP types = getListElement(sexpemission, "types");
                int currDim = length(currDims);
                Ds = new int[currDim];
                const char* typ = CHAR(STRING_ELT(types, currEmission));
                int currD;
                for (currD = 0; currD < LENGTH(currDims); currD++)
                {
                    // in R rooted as 1
                    Ds[currD] = INTEGER(VECTOR_ELT(currDims, currEmission))[0]-1 ;
                }

                SEXP t=R_NilValue;
                PROTECT(t = GET_SLOT(VECTOR_ELT(myInpedendentEmissions, currEmission), install("dim")));
                GET_SLOT(VECTOR_ELT(myInpedendentEmissions, currEmission), install("dim"));
                int dim =INTEGER(t)[0];
                UNPROTECT(1);

                EmissionFunction** myEmissionsB = NULL;
                myEmissionsB = RGETEMISSION(VECTOR_ELT(myInpedendentEmissions, currEmission), dim, sexpk,Ds, typ);

                combinedDistributions.push_back(myEmissionsB);

            }
            myEmissions = createJointlyIndependent(combinedDistributions, D, sexpk, T, nsample);
        }
        else
        {
            error("Unknown emission function specified: %s\n", type);
        }
        return myEmissions;

    }

    SEXP prepareEmission(const char* type, SEXP sexpfixedEmission, SEXP sexpemission, EmissionFunction** myEmissions, int K)
    {
        SEXP sexpemissionParam;
        const char* gauss = "Gaussian";
        const char* independent = "Independent";
        const char* jointlyindependent = "JointlyIndependent";

        if(LENGTH(sexpfixedEmission) == 0)
        {
            if(strcmp(type, gauss) == 0)
            {
                sexpemissionParam = RPREPAREGAUSSPAR(myEmissions, K, 1);
            }
            else if(strcmp(type, jointlyindependent) == 0)
            {
                sexpemissionParam = RPREPAREJOINTLYINDEPENDENTPAR(myEmissions, K,
                                    getListElement(sexpemission, "types"));
            }
            else
            {
                sexpemissionParam = RPREPAREMIXEDBERNOULLIGAUSSPAR(myEmissions, K,
                                    getListElement(sexpemission, "types"),
                                    getListElement(sexpemission, "order"));
            }
        }
        else
        {
            PROTECT(sexpemissionParam = NEW_LIST(0));
        }
        return sexpemissionParam;
    }

    SEXP RHMMVITERBI(SEXP sexpobs, SEXP sexppi, SEXP sexpA, SEXP sexpemission, SEXP sexptype, SEXP sexpdim, SEXP sexpk, SEXP sexpverbose, SEXP sexpfixedEmission)
    {
        // memory allocation
        if(DEBUG_MEMORY)
        {
            Rprintf("\n### Constructing C++ objects from R input ###\n\n");
        }
        int D = INTEGER(sexpdim)[0];
        int sumD, p;
        sumD = 0;
        for (p = 0; p < length(sexpdim); p++)
        {
            sumD = sumD + INTEGER(sexpdim)[p];
        }
        D = sumD;

        int i,j,k,n,t;
        int K = INTEGER(sexpk)[0];

        // parse observations
        int nsample = length(sexpobs);

        int* T = NULL;
        if(nsample > 0)
        {
            T = (int*)malloc(sizeof(int)*nsample);
            for(n=0; n<nsample; n++)
            {
                T[n] = INTEGER(getAttrib(coerceVector(VECTOR_ELT(sexpobs, n), REALSXP), R_DimSymbol))[0];
            }
        }

        if(DEBUG_MEMORY)
        {
            Rprintf("new->T; (%d bytes)\n", sizeof(int)*nsample);
        }
        double*** obs = RGETOBS(sexpobs, T, nsample, D);
        int** isNaN = whichNaN(obs, nsample, T, D);
        // start Vector for type Gaussian
        int* start_d = (int*)malloc(sizeof(int)*D);
        int o;
        for (o = 0; o < D; o++)
        {
            start_d[o] = o;
        }

        if(nsample == 0)
        {
            nsample = LENGTH(sexpfixedEmission);
            T = (int*)malloc(sizeof(int)*nsample);
            for(n=0; n<nsample; n++)
            {
                T[n] = INTEGER(getAttrib(coerceVector(VECTOR_ELT(sexpfixedEmission, n), REALSXP), R_DimSymbol))[0];
            }
        }

        InitialProbability* initProb = RGETINITPROB(sexppi, K);
        TransitionMatrix* transMat = RGETTRANSMAT(sexpA, K);

        EmissionFunction** myEmissions = NULL;
        const char* type = CHAR(STRING_ELT(sexptype,0));

        if(LENGTH(sexpfixedEmission) == 0)
        {
            myEmissions = getEmission(type, sexpemission, sexpk, start_d, nsample, T, K, D);
        }

        HMM* myHMM = createHMM(0, K, initProb, transMat, myEmissions);

        int verbose = INTEGER(sexpverbose)[0];
        int** S = (int**)malloc(sizeof(int*)*nsample);
        for(n=0; n<nsample; n++)
        {
            S[n] = (int*)malloc(sizeof(int)*T[n]);
        }

        double*** fixedEmission = NULL;
        if(LENGTH(sexpfixedEmission) > 0)
        {
            fixedEmission = (double***)malloc(sizeof(double**)*nsample);
            for(n=0; n<nsample; n++)
            {
                fixedEmission[n] = (double**)malloc(sizeof(double*)*K);
                for(i=0; i<K; i++)
                {
                    fixedEmission[n][i] = (double*)malloc(sizeof(double)*T[n]);
                    for(t=0; t<T[n]; t++)
                    {
                        fixedEmission[n][i][t] = REAL(coerceVector(VECTOR_ELT(sexpfixedEmission, n), REALSXP))[t+T[n]*i];
                    }
                }
            }
        }
        myHMM->Viterbi(S, obs, nsample, T, verbose, isNaN, fixedEmission);

        SEXP viterbi, tempvit;

        PROTECT(viterbi = NEW_LIST(nsample));
        for(n=0; n<nsample; n++)
        {
            PROTECT(tempvit = NEW_INTEGER(T[n]));
            for(i=0; i<T[n]; i++)
            {
                INTEGER_POINTER(tempvit)[i] = S[n][i]+1;
            }
            SET_ELEMENT(viterbi, n, tempvit);
            UNPROTECT(1);
        }

        delete myHMM;

        for(n=0; n<nsample; n++)
        {
            free(S[n]);
        }
        free(S);

        if(fixedEmission == NULL)
        {
            for(n=0; n<nsample; n++)
            {
                for(t=0; t<T[n]; t++)
                {
                    free(obs[n][t]);
                }
                free(obs[n]);
                free(isNaN[n]);
            }
            free(obs);
            free(isNaN);
        }

        free(T);

        UNPROTECT(1);

        return viterbi;
    }

    SEXP RHMMFit(SEXP sexpobs, SEXP sexppi, SEXP sexpA, SEXP sexpemission, SEXP sexptype, SEXP sexpdim, SEXP sexpregularize, SEXP sexpk, SEXP sexpmaxIters, SEXP sexpparallel, SEXP sexpflags, SEXP sexpstate2flag, SEXP sexpcouples, SEXP sexprevop, SEXP sexpverbose, SEXP sexpupdateTransMat, SEXP sexpfixedEmission, SEXP bidirOptimParams, SEXP emissionPrior, SEXP sexpeffectivezero, SEXP sepconvergence, SEXP sexpincrementalEM)
    {

        // memory allocation
        if(DEBUG_MEMORY)
        {
            Rprintf("\n### Constructing C++ objects from R input ###\n\n");
        }
        int D = INTEGER(sexpdim)[0];
        int sumD, p;
        sumD = 0;
        for (p = 0; p < length(sexpdim); p++)
        {
            sumD = sumD + INTEGER(sexpdim)[p];
        }
        D = sumD;
        double regularize = REAL(sexpregularize)[0];

        int i,j,k,n,t;
        int K = INTEGER(sexpk)[0];
        int ncores = INTEGER(sexpparallel)[0];
        if(RUNOPENMPVERSION == 0) {
            ncores = 1;
	}

        // parse observations
        int nsample = length(sexpobs);
        int incrementalEM = INTEGER(sexpincrementalEM)[0];

        int* T = NULL;
        if(nsample > 0)
        {
            T = (int*)malloc(sizeof(int)*nsample);
            for(n=0; n<nsample; n++)
            {
                T[n] = INTEGER(getAttrib(coerceVector(VECTOR_ELT(sexpobs, n), REALSXP), R_DimSymbol))[0];
            }
        }

        if(DEBUG_MEMORY)
        {
            Rprintf("new->T; (%d bytes)\n", sizeof(int)*nsample);
        }
        double*** obs = RGETOBS(sexpobs, T, nsample, D);
        int** isNaN = whichNaN(obs, nsample, T, D);
        // start Vector for type Gaussian
        int* start_d = (int*)malloc(sizeof(int)*D);
        int o;
        for (o = 0; o < D; o++)
        {
            start_d[o] = o;
        }

        if(nsample == 0)
        {
            nsample = LENGTH(sexpfixedEmission);
            T = (int*)malloc(sizeof(int)*nsample);
            for(n=0; n<nsample; n++)
            {
                T[n] = INTEGER(getAttrib(coerceVector(VECTOR_ELT(sexpfixedEmission, n), REALSXP), R_DimSymbol))[0];
            }
        }

        int *couples = NULL;
        RGETCOUPLES(sexpcouples, &couples, K);

        InitialProbability* initProb = RGETINITPROB(sexppi, K);
        TransitionMatrix* transMat = RGETTRANSMAT(sexpA, K);

        EmissionFunction** myEmissions = NULL;
        const char* type = CHAR(STRING_ELT(sexptype,0));

        if(LENGTH(sexpfixedEmission) == 0)
        {
            myEmissions = getEmission(type, sexpemission, sexpk, start_d, nsample, T, K, D);
        }

        // parallel
        HMM* myHMM = createHMM(0, K, initProb, transMat, myEmissions);

        int maxIters = INTEGER(sexpmaxIters)[0];

        if(DEBUG_MEMORY)
        {
            Rprintf("\n### Object construction complete. Getting ready to rumble... ###\n\n");
        }
        if(DEBUG)
        {
            Rprintf("\n### Fitting HMM model parameters ###\n\n");
        }

        int **flags = NULL;
        int *state2flag = NULL;
        RGETFLAGS(sexpflags, sexpstate2flag, &flags, &state2flag, nsample, T, K);

        int verbose = INTEGER(sexpverbose)[0];
        int* revop = NULL;

        if(LENGTH(sexprevop) > 0)
        {
            if(DEBUG)
            {
                Rprintf("reverse.observations = [");
            }
            revop = (int*)malloc(sizeof(int)*D);
            for(i=0; i<LENGTH(sexprevop); i++)
            {
                revop[i] = INTEGER(sexprevop)[i];
                if(DEBUG)
                {
                    Rprintf("%d ", revop[i]);
                }
            }
            if(DEBUG)
            {
                Rprintf("]\n");
            }
        }

        int updateTransMat = INTEGER(sexpupdateTransMat)[0];

        double*** fixedEmission = NULL;
        if(LENGTH(sexpfixedEmission) > 0)
        {
            fixedEmission = (double***)malloc(sizeof(double**)*nsample);
            for(n=0; n<nsample; n++)
            {
                fixedEmission[n] = (double**)malloc(sizeof(double*)*K);
                for(i=0; i<K; i++)
                {
                    //t+T[n]*d
                    fixedEmission[n][i] = (double*)malloc(sizeof(double)*T[n]);
                    for(t=0; t<T[n]; t++)
                    {
                        fixedEmission[n][i][t] = REAL(coerceVector(VECTOR_ELT(sexpfixedEmission, n), REALSXP))[t+T[n]*i];
                    }
                }
            }
        }

        double effective_zero = REAL(sexpeffectivezero)[0];
        double convergence = REAL(sepconvergence)[0];

        list<double> loglik = myHMM->BaumWelch(obs, T, nsample, maxIters, flags, state2flag, couples, revop, verbose, updateTransMat, isNaN, fixedEmission, bidirOptimParams, emissionPrior, ncores, effective_zero, convergence, incrementalEM);

        if(DEBUG)
        {
            Rprintf("\n\n### Parameter fitting complete. ###\n\n");
        }

        if(DEBUG)
        {
            Rprintf("\n### Creating R output...");
        }

        SEXP sexpemissionParam = prepareEmission(type, sexpfixedEmission, sexpemission, myEmissions, K);

        SEXP sexpinitProb =  RPREPAREPI(initProb);
        SEXP sexptransMat = RPREPARETRANSMAT(transMat);

        SEXP log_likelihood, HMMFIT, wnames;

        double curr_ll;
        PROTECT(log_likelihood = NEW_NUMERIC(loglik.size()));
        list<double>::iterator start, stop;
        start = loglik.begin();
        stop = loglik.end();
        i=0;
        for(list<double>::iterator it = start; it != stop; it++)
        {
            curr_ll = *it;
            NUMERIC_POINTER(log_likelihood)[i++] = curr_ll;
        }

        PROTECT(HMMFIT = NEW_LIST(4));
        PROTECT(wnames = NEW_CHARACTER(4));

        SET_STRING_ELT(wnames, 0, mkChar("loglik"));
        SET_STRING_ELT(wnames, 1, mkChar("initProb"));
        SET_STRING_ELT(wnames, 2, mkChar("transMat"));
        SET_STRING_ELT(wnames, 3, mkChar("emission"));
        SET_NAMES(HMMFIT, wnames);
        UNPROTECT(1);

        SET_ELEMENT(HMMFIT, 0, log_likelihood);
        SET_ELEMENT(HMMFIT, 1, sexpinitProb);
        SET_ELEMENT(HMMFIT, 2, sexptransMat);
        SET_ELEMENT(HMMFIT, 3, sexpemissionParam);
        UNPROTECT(2);

        if(LENGTH(sexpfixedEmission) > 0)
        {
            UNPROTECT(1);
        }

        if(DEBUG)
        {
            Rprintf("done. ###\n\n");
        }

        if(DEBUG_MEMORY)
        {
            Rprintf("\n### C++ memory deallocation ###\n\n");
        }
        delete myHMM;

        int mem = 0;

        if(fixedEmission == NULL)
        {
            for(n=0; n<nsample; n++)
            {
                free(isNaN[n]);
                for(t=0; t<T[n]; t++)
                {
                    free(obs[n][t]);
                    mem += sizeof(double)*D;
                }
                free(obs[n]);
                mem += sizeof(double*)*T[n];
            }
            free(obs);
            free(isNaN);
            mem += sizeof(double**)*nsample;
            if(DEBUG_MEMORY)
            {
                Rprintf("delete->observation matrix; (%d bytes)\n", mem);
            }
        }

        free(T);
        if(DEBUG_MEMORY)
        {
            Rprintf("delete->T (%d bytes);\n", sizeof(int)*nsample);
        }
        if(LENGTH(sexprevop) > 0)
        {
            free(revop);
        }
        if(couples != NULL)
        {
            free(couples);
        }

        if(DEBUG_MEMORY)
        {
            Rprintf("\n### deallocation complete. ###\n\n");
        }

        loglik.clear();
        return HMMFIT;

    }

    SEXP RGETPOSTERIOR(SEXP sexpobs, SEXP sexppi, SEXP sexpA, SEXP sexpemission, SEXP sexptype, SEXP sexpdim, SEXP sexpk, SEXP sexpverbose, SEXP sexpfixedEmission, SEXP sexpncores, SEXP sexpflags, SEXP sexpstate2flag)
    {

        if(DEBUG_MEMORY)
        {
            Rprintf("\n### Constructing C++ objects from R input ###\n\n");
        }

        // memory allocation
        int D = INTEGER(sexpdim)[0];
        int sumD, p;
        sumD = 0;
        for (p = 0; p < length(sexpdim); p++)
        {
            sumD = sumD + INTEGER(sexpdim)[p];
        }
        D = sumD;

        int i,j,k,n,t;
        int K = INTEGER(sexpk)[0];

        // parse observations
        int nsample = length(sexpobs);

        int* T = NULL;
        if(nsample > 0)
        {
            T = (int*)malloc(sizeof(int)*nsample);
            for(n=0; n<nsample; n++)
            {
                T[n] = INTEGER(getAttrib(coerceVector(VECTOR_ELT(sexpobs, n), REALSXP), R_DimSymbol))[0];
            }
        }

        if(DEBUG_MEMORY)
        {
            Rprintf("new->T; (%d bytes)\n", sizeof(int)*nsample);
        }
        double*** obs = RGETOBS(sexpobs, T, nsample, D);
        int** isNaN = whichNaN(obs, nsample, T, D);
        // start Vector for type Gaussian
        int* start_d = (int*)malloc(sizeof(int)*D);
        int o;
        for (o = 0; o < D; o++)
        {
            start_d[o] = o;
        }

        if(nsample == 0)
        {
            nsample = LENGTH(sexpfixedEmission);
            T = (int*)malloc(sizeof(int)*nsample);
            for(n=0; n<nsample; n++)
            {
                T[n] = INTEGER(getAttrib(coerceVector(VECTOR_ELT(sexpfixedEmission, n), REALSXP), R_DimSymbol))[0];
            }
        }
        int ncores = INTEGER(sexpncores)[0];

        InitialProbability* initProb = RGETINITPROB(sexppi, K);
        TransitionMatrix* transMat = RGETTRANSMAT(sexpA, K);

        EmissionFunction** myEmissions = NULL;
        const char* type = CHAR(STRING_ELT(sexptype,0));

        if(LENGTH(sexpfixedEmission) == 0)
        {
            myEmissions = getEmission(type, sexpemission, sexpk, start_d, nsample, T, K, D);
        }

        HMM* myHMM = createHMM(0, K, initProb, transMat, myEmissions);

        if(DEBUG_MEMORY)
        {
            Rprintf("\n### Object construction complete. Getting ready to rumble... ###\n\n");
        }
        if(DEBUG)
        {
            Rprintf("\n### Fitting HMM model parameters ###\n\n");
        }

        int verbose = INTEGER(sexpverbose)[0];

        double*** fixedEmission = NULL;
        if(LENGTH(sexpfixedEmission) > 0)
        {
            fixedEmission = (double***)malloc(sizeof(double**)*nsample);
            for(n=0; n<nsample; n++)
            {
                fixedEmission[n] = (double**)malloc(sizeof(double*)*K);
                for(i=0; i<K; i++)
                {
                    fixedEmission[n][i] = (double*)malloc(sizeof(double)*T[n]);
                    for(t=0; t<T[n]; t++)
                    {
                        fixedEmission[n][i][t] = REAL(coerceVector(VECTOR_ELT(sexpfixedEmission, n), REALSXP))[t+T[n]*i];
                    }
                }
            }
        }

        // memory allocation
        double** alpha = NULL;
        double** beta = NULL;
        double** gamma = NULL;
        double** emissionProb = NULL;
        double* c = NULL;
        int maxLen;
        maxLen = 0;
        for(n=0; n<nsample; n++)
        {
            if(maxLen < T[n])
            {
                maxLen = T[n];
            }
        }
        allocateMemAlpha(&alpha, maxLen, K);
        allocateMemBeta(&beta, maxLen, K);
        allocateMemRescFac(&c, maxLen, K);
        allocateMemGamma(&gamma, maxLen, K);
        allocateMemEmissionProb(&emissionProb, maxLen, K);

        int **flags = NULL;
        int *state2flag = NULL;
        RGETFLAGS(sexpflags, sexpstate2flag, &flags, &state2flag, nsample, T, K);

        SEXP sexpgamma, sexpcurrGAMMA;
        PROTECT(sexpgamma = NEW_LIST(nsample));
        for(n=0; n<nsample; n++)
        {
            if(fixedEmission == NULL)
            {
                myHMM->calcEmissionProbs(obs, emissionProb, T, n, flags, state2flag, NULL, isNaN, 1, verbose);
            }
            else
            {
                for(i=0; i<K; i++)
                {
                    for(t=0; t<T[n]; t++)
                    {
                        emissionProb[i][t] = fixedEmission[n][i][t];

                    }
                }
            }
            myHMM->getAlphaBeta(obs, alpha, beta, c, emissionProb, T, n, 1, 0, verbose);
            myHMM->getGamma(alpha, beta, c,  emissionProb, gamma, T, n, ncores, -1, verbose);
            PROTECT(sexpcurrGAMMA = NEW_NUMERIC(T[n]*K));
            for(i=0; i<K; i++)
            {
                for(t=0; t<T[n]; t++)
                {
                    //c[t];
                    NUMERIC_POINTER(sexpcurrGAMMA)[t+T[n]*i] = gamma[t][i];
                }
            }
            SET_ELEMENT(sexpgamma, n, sexpcurrGAMMA);
        }

        UNPROTECT(1+nsample);

        deallocateMemAlpha(alpha, maxLen, K);
        deallocateMemBeta(beta, maxLen, K);
        deallocateMemRescFac(c, maxLen, K);
        deallocateMemGamma(gamma, maxLen, K);
        deallocateMemEmissionProb(emissionProb, maxLen, K);

        delete myHMM;
        if(fixedEmission == NULL)
        {
            for(n=0; n<nsample; n++)
            {
                free(isNaN[n]);
                for(t=0; t<T[n]; t++)
                {
                    free(obs[n][t]);
                }
                free(obs[n]);

            }
            free(isNaN);
            free(obs);
        }

        if(fixedEmission != NULL)
        {
            for(n=0; n<nsample; n++)
            {
                for(i=0; i<K; i++)
                {
                    free(fixedEmission[n][i]);
                }
                free(fixedEmission[n]);
            }
            free(fixedEmission);
        }
        free(T);

        return sexpgamma;

    }

}
