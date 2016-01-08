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
        {"RGETLOGLIK", (DL_FUNC)&RGETLOGLIK, 12},
        {NULL, NULL, 0}
    };

    void R_init_STAN(DllInfo *info)
    {
        R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    }

    EmissionFunction** RGETMULTGAUSS(SEXP sexpmu, SEXP sexpsigma, int D, SEXP sexpk,int* start, int updateCov, int sharedCov)
    {
//int D = INTEGER(sexpdim)[0];
        int parallel = 0;

        int i,j,k;

        int K = INTEGER(sexpk)[0];

        EmissionFactory* factory = createEmissionFactory(MULTIVARIATEGAUSSIAN);
        EmissionFunction **HMMEmissionFunctions = allocateEmissionFunctionVector(K);

        if(DEBUG)
        {
            Rprintf("Getting mutlivariate Gaussian: \n");
            Rprintf("dim=%d \nstart=c(%d", D, start[0]);
            for(i=1; i<D; i++)
            {
                Rprintf(", %d", start[i]);
            }
            Rprintf(")\n");
        }

        for(k=0; k<K; k++)
        {
//	Rprintf("\tstate=%d: ", k);
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
            HMMEmissionFunctions[k] = factory->createEmissionFunction(new ParamContainerEmissions(mu, sigma, 0, D,start, updateCov, sharedCov), parallel);
        }

        delete factory;

        return HMMEmissionFunctions;
    }

    EmissionFunction** RGETBERNOULLI2(SEXP sexpbernoullip, int D, SEXP sexpk, int* start)
    {
        int i,j,k;
//	Rprintf("getting bernoulli\n");
        int K = INTEGER(sexpk)[0];

        EmissionFactory* factory = createEmissionFactory(BERNOULLI);
        EmissionFunction **HMMEmissionFunctions = allocateEmissionFunctionVector(K);
//Rprintf("Ber-D=%d\n", D);
        for(k=0; k<K; k++)
        {
//for p's
            double p1 = REAL(coerceVector(VECTOR_ELT(sexpbernoullip, k), REALSXP))[0];
// create new constructor in ParamContainerEmissions (p* and dim(int) )
            HMMEmissionFunctions[k] = factory->createEmissionFunction(new ParamContainerEmissions(p1, D, start), 0);
        }

        delete factory;

        return HMMEmissionFunctions;
    }

    EmissionFunction** RGETPOISSON(SEXP sexppoissonlambda, int D, SEXP sexpk, int* start)
    {
        int i,j,k;
        int K = INTEGER(sexpk)[0];

        EmissionFactory* factory = createEmissionFactory(POISSON);
        EmissionFunction **HMMEmissionFunctions = allocateEmissionFunctionVector(K);

        for(k=0; k<K; k++)
        {
//for p's
            double lambda = REAL(coerceVector(VECTOR_ELT(sexppoissonlambda, k), REALSXP))[0];
// create new constructor in ParamContainerEmissions (p* and dim(int) )
            HMMEmissionFunctions[k] = factory->createEmissionFunction(new ParamContainerEmissions(lambda, D, start, POISSON), 0);
        }

        delete factory;

        return HMMEmissionFunctions;
    }

    EmissionFunction** RGETMULTINOMIAL(SEXP sexpmultiP, SEXP sexpmultiRevCompl, int D, SEXP sexpk, int* start, int* state2flag)
    {
        int i,j,d,k;
        int K = INTEGER(sexpk)[0];

        EmissionFactory* factory = createEmissionFactory(MULTINOMIAL);
        EmissionFunction **HMMEmissionFunctions = allocateEmissionFunctionVector(K);
        for(k=0; k<K; k++)
        {
//for p's
            double* p = (double*)malloc(sizeof(double)*D);
            for(d=0; d<D; d++)
            {
                p[d] = REAL(coerceVector(VECTOR_ELT(sexpmultiP, k), REALSXP))[d];
            }
            int* revcomp = (int*)malloc(sizeof(int)*D);
            for(d=0; d<D; d++)
            {
                revcomp[d] = INTEGER(sexpmultiRevCompl)[d]-1;
//		printf("%d ", revcomp[d]);
            }
// create new constructor in ParamContainerEmissions (p* and dim(int) )
            if(state2flag != NULL)
            {
//	printf("state2flag[%d]=%d\n", k, state2flag[k]);

                HMMEmissionFunctions[k] = factory->createEmissionFunction(new ParamContainerEmissions(p, revcomp, D, start, state2flag[k]), 0);
            }
            else
            {
                HMMEmissionFunctions[k] = factory->createEmissionFunction(new ParamContainerEmissions(p, revcomp, D, start, -100), 0);
            }
        }

        delete factory;

        return HMMEmissionFunctions;
    }

    EmissionFunction** RGETNEGATIVEBINOMIAL(SEXP sexpMU, SEXP sexpSIZE, SEXP sexpSIZEFAC, SEXP sexpPI, int D, SEXP sexpk, int* start, double*** observations, int* T, int nsample, SEXP uniqueCountSplit, int* revop)
    {
        int i,j,d,k,n;
        int K = INTEGER(sexpk)[0];

        EmissionFactory* factory = createEmissionFactory(NEGATIVEBINOMIAL);
        EmissionFunction **HMMEmissionFunctions = allocateEmissionFunctionVector(K);
        for(k=0; k<K; k++)
        {
//for p's
            double mu_nb = REAL(coerceVector(VECTOR_ELT(sexpMU, k), REALSXP))[0];
            double mu_size = REAL(coerceVector(VECTOR_ELT(sexpSIZE, k), REALSXP))[0];

            int myNSamp = LENGTH(VECTOR_ELT(sexpSIZEFAC, k));
            double* sizeFactor_nb = (double*)malloc(sizeof(double)*myNSamp);
            for(i=0; i<myNSamp; i++)
            {
                sizeFactor_nb[i] = REAL(coerceVector(VECTOR_ELT(sexpSIZEFAC, k), REALSXP))[i];
//		printf("sF[%d]=%f\n", i, sizeFactor_nb[i]);
            }
//double* sizeFactor_nb = REAL(coerceVector(VECTOR_ELT(sexpSIZEFAC, k), REALSXP))[0];
            double pi_nb = REAL(coerceVector(VECTOR_ELT(sexpPI, k), REALSXP))[0];

// create new constructor in ParamContainerEmissions (p* and dim(int) )
            HMMEmissionFunctions[k] = factory->createEmissionFunction(new ParamContainerEmissions(mu_nb, mu_size, sizeFactor_nb, pi_nb, D, start,uniqueCountSplit), 0);
            if(observations != NULL)
            {
                HMMEmissionFunctions[k]->getParameter()->initUniqueObsProb(observations, nsample, T, revop);
                double** upobs = HMMEmissionFunctions[k]->getParameter()->getUniqueObsProb();
                int** ulens = HMMEmissionFunctions[k]->getParameter()->getUniqueLens();
                double* myval = (double*)malloc(sizeof(double)*1);

                for(n=0; n<nsample; n++)
                {
                    for(j=0; j<ulens[n][0]; j++)
                    {
                        if(upobs[n][j] != -1)
                        {
                            myval[0] = (double)j;
                            upobs[n][j] = HMMEmissionFunctions[k]->calcEmissionProbability(myval, -1, n);
                        }
                    }
                }

                free(myval);
            }
        }

        delete factory;

        return HMMEmissionFunctions;
    }

    EmissionFunction** RGETPOISSONLOGNORMAL(SEXP sexpMU, SEXP sexpSIGMA, SEXP sexpSIZEFAC, int D, SEXP sexpk, int* start, double*** observations, int* T, int nsample, SEXP uniqueCountSplit, int* revop)
    {
        int i,j,d,k,n;
        int K = INTEGER(sexpk)[0];

        EmissionFactory* factory = createEmissionFactory(POISSONLOGNORMAL);
        EmissionFunction **HMMEmissionFunctions = allocateEmissionFunctionVector(K);
        for(k=0; k<K; k++)
        {
//for p's
            double mu_pln = REAL(coerceVector(VECTOR_ELT(sexpMU, k), REALSXP))[0];
            double sigma_pln = REAL(coerceVector(VECTOR_ELT(sexpSIGMA, k), REALSXP))[0];
            double* sigma_pln_vec = NULL;
            int myNSamp = LENGTH(VECTOR_ELT(sexpSIZEFAC, k));
/*if(LENGTH(coerceVector(VECTOR_ELT(sexpSIGMA, k), REALSXP)) > 1) {
    sigma_pln_vec = (double*)malloc(sizeof(double)*myNSamp);
    for(i=0; i<myNSamp; i++) {
        sigma_pln_vec[i] = REAL(coerceVector(VECTOR_ELT(sexpSIGMA, k), REALSXP))[i];
    }

}
else {
    sigma_pln = REAL(coerceVector(VECTOR_ELT(sexpSIGMA, k), REALSXP))[0];
//	printf("sif=%f\n", sigma_pln);
}*/
            double* sizeFactor_pln = (double*)malloc(sizeof(double)*myNSamp);
            for(i=0; i<myNSamp; i++)
            {
                sizeFactor_pln[i] = REAL(coerceVector(VECTOR_ELT(sexpSIZEFAC, k), REALSXP))[i];
//printf("sF[%d]=%f\n", i, sizeFactor_pln[i]);
            }

// create new constructor in ParamContainerEmissions (p* and dim(int) )
            HMMEmissionFunctions[k] = factory->createEmissionFunction(new ParamContainerEmissions(mu_pln, sigma_pln, sigma_pln_vec, sizeFactor_pln, D, start, uniqueCountSplit), 0);
            if(observations != NULL)
            {
                HMMEmissionFunctions[k]->getParameter()->initUniqueObsProb(observations, nsample, T, revop);
                double** upobs = HMMEmissionFunctions[k]->getParameter()->getUniqueObsProb();
                int** ulens = HMMEmissionFunctions[k]->getParameter()->getUniqueLens();
                double* myval = (double*)malloc(sizeof(double)*1);

                for(n=0; n<nsample; n++)
                {
                    for(j=0; j<ulens[n][0]; j++)
                    {
                        if(upobs[n][j] != -1)
                        {
                            myval[0] = (double)j;
                            upobs[n][j] = HMMEmissionFunctions[k]->calcEmissionProbability(myval, -1, n);
                        }
                    }
                }

                free(myval);
            }
        }

        delete factory;

        return HMMEmissionFunctions;
    }

    EmissionFunction** RGETEMISSION(SEXP sexpparameters, int D, SEXP sexpk, int* start, const char* type, double*** observations, int* T, int nsample, SEXP uniqueCountSplit, int* revop, int* state2flag, int* couples)
    {
        EmissionFunction **HMMEmissionFunctions;
        const char* gaussian = "Gaussian";
        const char* bernoulli = "Bernoulli";
        const char* poisson = "Poisson";
        const char* multinomial = "Multinomial";
        const char* negativebinomial = "NegativeBinomial";
        const char* poissonlognormal = "PoissonLogNormal";

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
            HMMEmissionFunctions = RGETMULTGAUSS(getListElement(t, "mu"), getListElement(t, "cov"), D, sexpk, start, INTEGER(getListElement(t, "updateCov"))[0], INTEGER(getListElement(t, "sharedCov"))[0]);
        }
        else if(strcmp(type,poisson) == 0)
        {
            HMMEmissionFunctions = RGETPOISSON(getListElement(t, "lambda"), D, sexpk, start);
        }
        else if(strcmp(type,multinomial) == 0)
        {
            HMMEmissionFunctions = RGETMULTINOMIAL(getListElement(t, "p"), getListElement(t, "reverseComplementary"), D, sexpk, start, state2flag);
        }
        else if(strcmp(type,negativebinomial) == 0)
        {
            HMMEmissionFunctions = RGETNEGATIVEBINOMIAL(getListElement(t, "mu"), getListElement(t, "size"), getListElement(t, "sizeFactor"), getListElement(t, "pi"), D, sexpk, start, observations, T, nsample, getListElement(uniqueCountSplit, "NegativeBinomial"), revop);
        }
        else if(strcmp(type,poissonlognormal) == 0)
        {
            HMMEmissionFunctions = RGETPOISSONLOGNORMAL(getListElement(t, "mu"), getListElement(t, "sigma"), getListElement(t, "sizeFactor"), D, sexpk, start, observations, T, nsample, getListElement(uniqueCountSplit, "PoissonLogNormal"), revop);
        }

        UNPROTECT(1);

        return HMMEmissionFunctions;
    }

/*
 * BERNOULLI
 */

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

//Rprintf("Creating output for multivariate gaussian emission density.\n");
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
//Rprintf("%f ", myEmissions[k]->getParameter()->getGaussianSIGMA()[i][j]);
                }
//Rprintf("\n");
            }
            SET_ELEMENT(sigmaFit, k, currSIGMA);
//	Rprintf("\n");
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
            SET_STRING_ELT(wnames, 0, mkChar("mu"));
            SET_STRING_ELT(wnames, 1, mkChar("cov"));
            SET_STRING_ELT(wnames, 2, mkChar("invsigma"));
            SET_NAMES(emissionParam, wnames);
            UNPROTECT(1);
        }

        SET_ELEMENT(emissionParam, 0, muFit);
        SET_ELEMENT(emissionParam, 1, sigmaFit);  //sigmaFit);
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
//Rprintf("myEmission BernoulliP, %f ", myEmissions[k]->getParameter()->getBernoulliP());

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
//Rprintf("myEmission BernoulliP, %f \n", myEmissions[k]->getParameter()->getBernoulliP());
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

    SEXP RPREPAREPOISSONPAR(EmissionFunction** myEmissions, int K, int useNames)
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
                NUMERIC_POINTER(currState)[d] = myEmissions[k]->getParameter()->getPoissonLambda();
//Rprintf("myEmission BernoulliP, %f \n", myEmissions[k]->getParameter()->getBernoulliP());
            }
            SET_ELEMENT(p, k, currState);
        }

        if(useNames)
        {
            PROTECT(wname = NEW_CHARACTER(1));
            SET_STRING_ELT(wname, 0, mkChar("lambda"));
            SET_NAMES(p, wname);
            UNPROTECT(1);
        }

        UNPROTECT(K+1);

        return p;

    }

    SEXP RPREPAREMULTINOMIALPAR(EmissionFunction** myEmissions, int K, int useNames)
    {
        int k,d,j;
        int D = myEmissions[0]->getParameter()->getD();

        SEXP p, emissionParam, myComp, wname, currP;
        PROTECT(emissionParam = NEW_LIST(2));

        PROTECT(p = NEW_LIST(K));
        for(k=0; k<K; k++)
        {
            PROTECT(currP = NEW_NUMERIC(D));
            for(d=0; d<D; d++)
            {
                NUMERIC_POINTER(currP)[d] = myEmissions[k]->getParameter()->getMultinomialP()[d];

            }
            SET_ELEMENT(p, k, currP);
        }

        PROTECT(myComp = NEW_INTEGER(D));
        for(d=0; d<D; d++)
        {
//	printf("d=%\n", myEmissions[k]->getParameter()->getReverseComplementary()[d]+1);
            INTEGER_POINTER(myComp)[d] = myEmissions[0]->getParameter()->getReverseComplementary()[d]+1;

        }

        SET_ELEMENT(emissionParam, 0, p);
        SET_ELEMENT(emissionParam, 1, myComp);

        if(useNames)
        {
            PROTECT(wname = NEW_CHARACTER(1));
            SET_STRING_ELT(wname, 0, mkChar("p"));
            SET_STRING_ELT(wname, 1, mkChar("reverseComplementary"));
            SET_NAMES(emissionParam, wname);
            UNPROTECT(1);
        }

        UNPROTECT(K+3);

        return emissionParam;

    }

    SEXP RPREPARENEGATIVEBINOMIALPAR(EmissionFunction** myEmissions, int K, int useNames)
    {
        int k,d,j;
        int D = myEmissions[0]->getParameter()->getD();

        SEXP mu_nb, size_nb, sizeFactor, pi, emissionParam, wname, currMU, currSize, currSizeFactor, currPi;
        PROTECT(emissionParam = NEW_LIST(4));

        PROTECT(mu_nb = NEW_LIST(K));
        PROTECT(size_nb = NEW_LIST(K));
        PROTECT(sizeFactor = NEW_LIST(K));
        PROTECT(pi = NEW_LIST(K));

        for(k=0; k<K; k++)
        {
            PROTECT(currMU = NEW_NUMERIC(D));
            PROTECT(currSize = NEW_NUMERIC(D));
            PROTECT(currSizeFactor = NEW_NUMERIC(D));
            PROTECT(currPi = NEW_NUMERIC(D));

            for(d=0; d<D; d++)
            {
//Rprintf("d=%d\n",d);
                NUMERIC_POINTER(currMU)[d] = myEmissions[k]->getParameter()->getMuNB();
                NUMERIC_POINTER(currSize)[d] = myEmissions[k]->getParameter()->getSizeNB();
                NUMERIC_POINTER(currSizeFactor)[d] = myEmissions[k]->getParameter()->getSizeFactorNB()[0];
                NUMERIC_POINTER(currPi)[d] = myEmissions[k]->getParameter()->getPiNB();

            }
            SET_ELEMENT(mu_nb, k, currMU);
            SET_ELEMENT(size_nb, k, currSize);
            SET_ELEMENT(sizeFactor, k, currSizeFactor);
            SET_ELEMENT(pi, k, currPi);
        }

        SET_ELEMENT(emissionParam, 0, mu_nb);
        SET_ELEMENT(emissionParam, 1, size_nb);
        SET_ELEMENT(emissionParam, 2, sizeFactor);
        SET_ELEMENT(emissionParam, 3, pi);

        if(useNames)
        {
            PROTECT(wname = NEW_CHARACTER(4));
            SET_STRING_ELT(wname, 0, mkChar("mu"));
            SET_STRING_ELT(wname, 1, mkChar("size"));
            SET_STRING_ELT(wname, 2, mkChar("sizeFactor"));
            SET_STRING_ELT(wname, 3, mkChar("pi"));
            SET_NAMES(emissionParam, wname);
            UNPROTECT(2);
        }

        UNPROTECT(2*K+7);

        return emissionParam;

    }

    SEXP RPREPAREPOISSONLOGNORMALPAR(EmissionFunction** myEmissions, int K, int useNames)
    {
        int k,d,j;
        int D = myEmissions[0]->getParameter()->getD();

        SEXP mu_pln, sigma_pln, emissionParam, wname, currMU, currSigma, sizeFactor, currSizeFactor;
        PROTECT(emissionParam = NEW_LIST(3));

        PROTECT(mu_pln = NEW_LIST(K));
        PROTECT(sigma_pln = NEW_LIST(K));
        PROTECT(sizeFactor = NEW_LIST(K));

        for(k=0; k<K; k++)
        {
            PROTECT(currMU = NEW_NUMERIC(D));
            PROTECT(currSigma = NEW_NUMERIC(D));
            PROTECT(currSizeFactor = NEW_NUMERIC(D));

            for(d=0; d<D; d++)
            {
//Rprintf("d=%d\n",d);
                NUMERIC_POINTER(currMU)[d] = myEmissions[k]->getParameter()->getMuPoiLog();
                NUMERIC_POINTER(currSigma)[d] = myEmissions[k]->getParameter()->getSigmaPoiLog();
                NUMERIC_POINTER(currSizeFactor)[d] = myEmissions[k]->getParameter()->getSizeFactorPoiLog()[0];

            }
            SET_ELEMENT(mu_pln, k, currMU);
            SET_ELEMENT(sigma_pln, k, currSigma);
            SET_ELEMENT(sizeFactor, k, currSizeFactor);

        }

        SET_ELEMENT(emissionParam, 0, mu_pln);
        SET_ELEMENT(emissionParam, 1, sigma_pln);
        SET_ELEMENT(emissionParam, 2, sizeFactor);

        if(useNames)
        {
            PROTECT(wname = NEW_CHARACTER(4));
            SET_STRING_ELT(wname, 0, mkChar("mu"));
            SET_STRING_ELT(wname, 1, mkChar("sigma"));
            SET_STRING_ELT(wname, 1, mkChar("sizeFactor"));
            SET_NAMES(emissionParam, wname);
            UNPROTECT(1);
        }

        UNPROTECT(3*K+4);

        return emissionParam;

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
        else if(strcmp(type,"Poisson") == 0)
        {
            return RPREPAREPOISSONPAR(myEmissions, K, useNames);
        }
        else if(strcmp(type,"Multinomial") == 0)
        {
            return RPREPAREMULTINOMIALPAR(myEmissions, K, useNames);
        }
        else if(strcmp(type,"NegativeBinomial") == 0)
        {
            return RPREPARENEGATIVEBINOMIALPAR(myEmissions, K, useNames);
        }
        else if(strcmp(type,"PoissonLogNormal") == 0)
        {
            return RPREPAREPOISSONLOGNORMALPAR(myEmissions, K, useNames);
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
                free(currEmission);
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

        listEF.clear();
        UNPROTECT(3+LENGTH(types));

        return output;                            //emissionList;
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
            int d = 0;
            for (pos = bernoulligauss.begin(); pos!=bernoulligauss.end(); pos++)
            {
////Rprintf("\nthis start: %d  dim %d\n\n", (**pos)->getParameter()->getStart(), (**pos)->getParameter()->getD());
                EmissionFunction* f = (*pos)[k];
                f->getParameter()->setCurrState(k);
                combinedPerState.push_back(f);
//allEmissions[d++] = f;
//	delete f;
            }
/*std::list<EmissionFunction*>::iterator pos1;
for (pos1 = combinedPerState.begin(); pos1!=combinedPerState.end(); pos1++) {
    ////Rprintf("\nthis start: %d  dim %d\n\n", (**pos)->getParameter()->getStart(), (**pos)->getParameter()->getD());
    EmissionFunction* f = (*pos1);
    delete f;
}*//*
    for(d=0; d<D; d++) {
        delete allEmissions[d];
    }*/
            HMMEmissionFunctions[k] = factory->createEmissionFunctionMixed(combinedPerState, new ParamContainerEmissions(D));
            HMMEmissionFunctions[k]->getParameter()->setDataVars(nsample, T);
            HMMEmissionFunctions[k]->getParameter()->setCurrState(k);

            list<EmissionFunction*>::iterator subEmission;
            std::list<EmissionFunction*> listEF = HMMEmissionFunctions[k]->getEmissionFunctionList();
// set gamma aux of "sub-emissions" to gamma aux of wrapper emission
            for (subEmission = listEF.begin(); subEmission != listEF.end(); subEmission++)
            {
                double **myGammaAux = HMMEmissionFunctions[k]->getParameter()->getGammaAux();
                ((*subEmission)->getParameter())->setDataVars(myGammaAux, nsample, T);
//	((*subEmission)->getParameter())->setDataVars(nsample, T);
            }
        }

        delete factory;
        return HMMEmissionFunctions;
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
//	Rprintf("%f ", observations[n][t][d]);
                    }
//Rprintf("\n");
                }
//Rprintf("\n\n");
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

    EmissionFunction** getEmission(const char* type, SEXP sexpemission, SEXP sexpk, int* start_d, int nsample, int* T, int K, int D, double*** observations, int* revop, int* state2flag, int* couples)
    {
        int i;
        EmissionFunction** myEmissions = NULL;
        const char* gauss = "Gaussian";
        const char* independent = "Independent";
        const char* jointlyindependent = "JointlyIndependent";
        const char* multinomial = "Multinomial";
        const char* negativebinomial = "NegativeBinomial";

        if(strcmp(type, gauss) == 0)
        {
                                                  //parallel
            myEmissions = RGETMULTGAUSS(getListElement(sexpemission, "mu"), getListElement(sexpemission, "cov"), D, sexpk, start_d, INTEGER(getListElement(sexpemission, "updateCov"))[0], INTEGER(getListElement(sexpemission, "sharedCov"))[0]);
            for(i=0; i<K; i++)
            {
                myEmissions[i]->getParameter()->setDataVars(nsample, T);
            }
        }
        else if(strcmp(type, multinomial) == 0)
        {
            myEmissions = RGETMULTINOMIAL(getListElement(sexpemission, "p"), getListElement(sexpemission, "reverseComplementary"), D, sexpk, start_d, state2flag);
            for(i=0; i<K; i++)
            {
                myEmissions[i]->getParameter()->setDataVars(nsample, T);
            }
        }
        else if(strcmp(type, jointlyindependent) == 0)
        {
                                                  //length(getListElement(sexpemission, "emissionDim"));
            int nemissions = LENGTH(getListElement(sexpemission, "emissions"));
            list<EmissionFunction**> combinedDistributions;
            SEXP uniqueCountSplit = getListElement(sexpemission, "mySplit");

            SEXP myInpedendentEmissions = getListElement(sexpemission, "emissions");
            int currEmission;
            for (currEmission = 0; currEmission < nemissions; currEmission++)
            {
                int* Ds = NULL;
                int ordersize;
                SEXP currDims = getListElement(sexpemission, "emissionDim");
                SEXP types = getListElement(sexpemission, "types");
                int currDim = length(currDims);
                Ds = new int[LENGTH(VECTOR_ELT(currDims, currEmission))];
                const char* typ = CHAR(STRING_ELT(types, currEmission));
                int currD;
//Rprintf("L=%d\n", LENGTH(currDims));
                for (currD = 0; currD < LENGTH(VECTOR_ELT(currDims, currEmission)); currD++)
                {
                                                  // in R rooted as 1
                    Ds[currD] = INTEGER(VECTOR_ELT(currDims, currEmission))[currD]-1 ;
//Rprintf(" Ds[%d]: %d \n ", currD, Ds[currD]);
                }

                SEXP t=R_NilValue;
                PROTECT(t = GET_SLOT(VECTOR_ELT(myInpedendentEmissions, currEmission), install("dim")));
                GET_SLOT(VECTOR_ELT(myInpedendentEmissions, currEmission), install("dim"));
                int dim =INTEGER(t)[0];
                UNPROTECT(1);

                EmissionFunction** myEmissionsB = NULL;

                myEmissionsB = RGETEMISSION(VECTOR_ELT(myInpedendentEmissions, currEmission), dim, sexpk,Ds, typ, observations, T, nsample, uniqueCountSplit, revop, state2flag, couples);
//	delete myEmissionsB[0];
                combinedDistributions.push_back(myEmissionsB);

            }
                                                  //parallel
            myEmissions = createJointlyIndependent(combinedDistributions, D, sexpk, T, nsample);
            combinedDistributions.clear();

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
        const char* multinomial = "Multinomial";
        const char* negativebinomial = "NegativeBinomial";
        const char* poissonlognormal = "PoissonLogNormal";

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
            else if(strcmp(type, multinomial) == 0)
            {
                sexpemissionParam = RPREPAREMULTINOMIALPAR(myEmissions, K, 1);
            }
            else if(strcmp(type, negativebinomial) == 0)
            {
                RPREPARENEGATIVEBINOMIALPAR(myEmissions, K, 1);
            }
            else if(strcmp(type, negativebinomial) == 0)
            {
                RPREPAREPOISSONLOGNORMALPAR(myEmissions, K, 1);
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
//Rprintf("T[n]=%d\n", T[n]);
            }
        }

        if(DEBUG_MEMORY)
        {
            Rprintf("new->T; (%d bytes)\n", sizeof(int)*nsample);
        }
        double*** obs = RGETOBS(sexpobs, T, nsample, D);
        int** isNaN = whichNaN(obs, nsample, T, D);
//start Vector for type Gaussian
        int* start_d = (int*)malloc(sizeof(int)*D);
        int o;
        for (o = 0; o < D; o++)
        {
            start_d[o] = o;
//Rprintf("init %d ", start_d[o]);
        }

        if(nsample == 0)
        {
            nsample = LENGTH(sexpfixedEmission);
            T = (int*)malloc(sizeof(int)*nsample);
            for(n=0; n<nsample; n++)
            {
                T[n] = INTEGER(getAttrib(coerceVector(VECTOR_ELT(sexpfixedEmission, n), REALSXP), R_DimSymbol))[0];
//Rprintf("T[n]=%d\n", T[n]);
            }
        }

        InitialProbability* initProb = RGETINITPROB(sexppi, K);
        TransitionMatrix* transMat = RGETTRANSMAT(sexpA, K);

        EmissionFunction** myEmissions = NULL;
        const char* type = CHAR(STRING_ELT(sexptype,0));

        if(LENGTH(sexpfixedEmission) == 0)
        {
            myEmissions = getEmission(type, sexpemission, sexpk, start_d, nsample, T, K, D, obs, NULL, NULL, NULL);
        }

                                                  //parallel
        HMM* myHMM = createHMM(0, K, initProb, transMat, myEmissions);

        int verbose = INTEGER(sexpverbose)[0];
        int** S = (int**)malloc(sizeof(int*)*nsample);
        for(n=0; n<nsample; n++)
        {
            S[n] = (int*)malloc(sizeof(int)*T[n]);
        }

/*double*** fixedEmission = NULL;
if(LENGTH(sexpfixedEmission) > 0) {
    fixedEmission = (double***)malloc(sizeof(double**)*nsample);
    for(n=0; n<nsample; n++) {
        fixedEmission[n] = (double**)malloc(sizeof(double*)*K);
        for(i=0; i<K; i++) {
            fixedEmission[n][i] = (double*)malloc(sizeof(double)*T[n]);
            for(t=0; t<T[n]; t++) {
                fixedEmission[n][i][t] = REAL(coerceVector(VECTOR_ELT(sexpfixedEmission, n), REALSXP))[i+t*K];
                if(t<3) {
                                        Rprintf("i=%d, t=%d, %f \n", i, t, fixedEmission[n][i][t]);}
}
if(t<3){
Rprintf("\n");}
}
}
}*/
        double*** fixedEmission = NULL;
        if(LENGTH(sexpfixedEmission) > 0)
        {
            fixedEmission = (double***)malloc(sizeof(double**)*nsample);
            for(n=0; n<nsample; n++)
            {
//	Rprintf("T=%d\n", T[n]);
                fixedEmission[n] = (double**)malloc(sizeof(double*)*K);
                for(i=0; i<K; i++)
                {
                                                  //t+T[n]*d
                    fixedEmission[n][i] = (double*)malloc(sizeof(double)*T[n]);
                    for(t=0; t<T[n]; t++)
                    {
                        fixedEmission[n][i][t] = REAL(coerceVector(VECTOR_ELT(sexpfixedEmission, n), REALSXP))[t+T[n]*i];
//	if(t<3) {
//		Rprintf("i=%d, t=%d, %f \n", i, t, fixedEmission[n][i][t]);//}
                    }
//if(t<3){
//Rprintf("\n");//}
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
        free(start_d);

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

// parse observations
        int nsample = length(sexpobs);
        int incrementalEM = INTEGER(sexpincrementalEM)[0];
//	Rprintf("%d\n", incrementalEM);

        int* T = NULL;
        if(nsample > 0)
        {
            T = (int*)malloc(sizeof(int)*nsample);
            for(n=0; n<nsample; n++)
            {
                T[n] = INTEGER(getAttrib(coerceVector(VECTOR_ELT(sexpobs, n), REALSXP), R_DimSymbol))[0];
//Rprintf("T[n]=%d\n", T[n]);
            }
        }

        if(DEBUG_MEMORY)
        {
            Rprintf("new->T; (%d bytes)\n", sizeof(int)*nsample);
        }
        double*** obs = RGETOBS(sexpobs, T, nsample, D);
        int** isNaN = whichNaN(obs, nsample, T, D);
//start Vector for type Gaussian
        int* start_d = (int*)malloc(sizeof(int)*D);
        int o;
        for (o = 0; o < D; o++)
        {
            start_d[o] = o;
//	Rprintf("init %d ", start_d[o]);
        }

        if(nsample == 0)
        {
            nsample = LENGTH(sexpfixedEmission);
            T = (int*)malloc(sizeof(int)*nsample);
            for(n=0; n<nsample; n++)
            {
                T[n] = INTEGER(getAttrib(coerceVector(VECTOR_ELT(sexpfixedEmission, n), REALSXP), R_DimSymbol))[0];
//Rprintf("T[n]=%d\n", T[n]);
            }
        }

        int *couples = NULL;
        RGETCOUPLES(sexpcouples, &couples, K);

        InitialProbability* initProb = RGETINITPROB(sexppi, K);
        TransitionMatrix* transMat = RGETTRANSMAT(sexpA, K);

        EmissionFunction** myEmissions = NULL;
        const char* type = CHAR(STRING_ELT(sexptype,0));

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

        int **flags = NULL;
        int *state2flag = NULL;
        RGETFLAGS(sexpflags, sexpstate2flag, &flags, &state2flag, nsample, T, K);

        if(LENGTH(sexpfixedEmission) == 0)
        {
            myEmissions = getEmission(type, sexpemission, sexpk, start_d, nsample, T, K, D, obs, revop, state2flag, couples);
        }

                                                  //parallel
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

        int verbose = INTEGER(sexpverbose)[0];

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
//if(t<3) {
//	Rprintf("i=%d, t=%d, %f \n", i, t, fixedEmission[n][i][t]);}
                    }
//if(t<3){
//Rprintf("\n");}
                }
            }
        }

        double effective_zero = REAL(sexpeffectivezero)[0];
        double convergence = REAL(sepconvergence)[0];

        list<double> loglik = myHMM->BaumWelch(obs, T, nsample, maxIters, flags, state2flag, couples, revop, verbose, updateTransMat, isNaN, fixedEmission, bidirOptimParams, emissionPrior, ncores, effective_zero, convergence, incrementalEM);

        double* dirScore = NULL;

        SEXP sexp_dirscore;
        if(couples != NULL)
        {
            dirScore = (double*)malloc(sizeof(double)*K);
            myHMM->getDirScore(dirScore, flags, state2flag, couples, revop, isNaN, obs, fixedEmission, K, nsample, T, ncores, effective_zero);
            PROTECT(sexp_dirscore = NEW_NUMERIC(K));
            for(i=0; i<K; i++)
            {
                NUMERIC_POINTER(sexp_dirscore)[i] = dirScore[i];
            }
        }
        else
        {
            PROTECT(sexp_dirscore = NEW_NUMERIC(0));
        }

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
//Rprintf("%f\n", *it);
            NUMERIC_POINTER(log_likelihood)[i++] = curr_ll;
        }

        PROTECT(HMMFIT = NEW_LIST(5));
        PROTECT(wnames = NEW_CHARACTER(5));

        SET_STRING_ELT(wnames, 0, mkChar("loglik"));
        SET_STRING_ELT(wnames, 1, mkChar("initProb"));
        SET_STRING_ELT(wnames, 2, mkChar("transMat"));
        SET_STRING_ELT(wnames, 3, mkChar("emission"));
        SET_STRING_ELT(wnames, 4, mkChar("dirScore"));
        SET_NAMES(HMMFIT, wnames);
        UNPROTECT(1);

        SET_ELEMENT(HMMFIT, 0, log_likelihood);
        SET_ELEMENT(HMMFIT, 1, sexpinitProb);
        SET_ELEMENT(HMMFIT, 2, sexptransMat);
        SET_ELEMENT(HMMFIT, 3, sexpemissionParam);
        SET_ELEMENT(HMMFIT, 4, sexp_dirscore);
        UNPROTECT(3);

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

        if(couples != NULL)
        {
            free(dirScore);
        }
//RFREEFLAGS(sexpflags, sexpstate2flag, flags, state2flag, nsample);
        loglik.clear();
        free(start_d);
        free(state2flag);

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
//start Vector for type Gaussian
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
//if(t<3) {
//	Rprintf("i=%d, t=%d, %f \n", i, t, fixedEmission[n][i][t]);}
                    }
//if(t<3){
//Rprintf("\n");}
                }
            }
        }

        if(nsample == 0)
        {
            nsample = LENGTH(sexpfixedEmission);
            T = (int*)malloc(sizeof(int)*nsample);
            for(n=0; n<nsample; n++)
            {
                T[n] = INTEGER(getAttrib(coerceVector(VECTOR_ELT(sexpfixedEmission, n), REALSXP), R_DimSymbol))[0];
//Rprintf("T[n]=%d\n", T[n]);
            }
        }
        int ncores = INTEGER(sexpncores)[0];

        InitialProbability* initProb = RGETINITPROB(sexppi, K);
        TransitionMatrix* transMat = RGETTRANSMAT(sexpA, K);
//		delete initProb;
//			delete transMat;

        EmissionFunction** myEmissions = NULL;
        const char* type = CHAR(STRING_ELT(sexptype,0));
        int* start_d = (int*)malloc(sizeof(int)*D);
        int o;
        for (o = 0; o < D; o++)
        {
            start_d[o] = o;
//Rprintf("init %d ", start_d[o]);
        }
        if(LENGTH(sexpfixedEmission) == 0)
        {
            myEmissions = getEmission(type, sexpemission, sexpk, start_d, nsample, T, K, D, obs, NULL, NULL, NULL);
        }
                                                  //parallel
        HMM* myHMM = createHMM(0, K, initProb, transMat, myEmissions);
//	delete myHMM;
        if(DEBUG_MEMORY)
        {
            Rprintf("\n### Object construction complete. Getting ready to rumble... ###\n\n");
        }
        if(DEBUG)
        {
            Rprintf("\n### Fitting HMM model parameters ###\n\n");
        }

        int verbose = INTEGER(sexpverbose)[0];

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

        int* couples = (int*)malloc(sizeof(int)*K);
        for(i=0; i<K; i++)
        {
            couples[i] = i;
        }

        SEXP sexpgamma, sexpcurrGAMMA;
        PROTECT(sexpgamma = NEW_LIST(nsample));
        for(n=0; n<nsample; n++)
        {
            if(fixedEmission == NULL)
            {
                myHMM->calcEmissionProbs(obs, emissionProb, T, n, flags, state2flag, NULL, isNaN, 1, verbose, couples);
//Rprintf("%f\n", emissionProb[0][0]);
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
//Rprintf("gamma done.\n");
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
        free(start_d);
        free(couples);

//	SEXP sexpgamma;
//	PROTECT(sexpgamma = NEW_LIST(0));
//	UNPROTECT(1);
        return sexpgamma;

    }

    SEXP RGETLOGLIK(SEXP sexpobs, SEXP sexppi, SEXP sexpA, SEXP sexpemission, SEXP sexptype, SEXP sexpdim, SEXP sexpk, SEXP sexpverbose, SEXP sexpfixedEmission, SEXP sexpncores, SEXP sexpflags, SEXP sexpstate2flag)
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
//start Vector for type Gaussian
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
//if(t<3) {
//	Rprintf("i=%d, t=%d, %f \n", i, t, fixedEmission[n][i][t]);}
                    }
//if(t<3){
//Rprintf("\n");}
                }
            }
        }

        if(nsample == 0)
        {
            nsample = LENGTH(sexpfixedEmission);
            T = (int*)malloc(sizeof(int)*nsample);
            for(n=0; n<nsample; n++)
            {
                T[n] = INTEGER(getAttrib(coerceVector(VECTOR_ELT(sexpfixedEmission, n), REALSXP), R_DimSymbol))[0];
//Rprintf("T[n]=%d\n", T[n]);
            }
        }
        int ncores = INTEGER(sexpncores)[0];

        InitialProbability* initProb = RGETINITPROB(sexppi, K);
        TransitionMatrix* transMat = RGETTRANSMAT(sexpA, K);
//		delete initProb;
//			delete transMat;

        EmissionFunction** myEmissions = NULL;
        const char* type = CHAR(STRING_ELT(sexptype,0));
        int* start_d = (int*)malloc(sizeof(int)*D);
        int o;
        for (o = 0; o < D; o++)
        {
            start_d[o] = o;
//Rprintf("init %d ", start_d[o]);
        }
        if(LENGTH(sexpfixedEmission) == 0)
        {
            myEmissions = getEmission(type, sexpemission, sexpk, start_d, nsample, T, K, D, obs, NULL, NULL, NULL);
        }
                                                  //parallel
        HMM* myHMM = createHMM(0, K, initProb, transMat, myEmissions);
//	delete myHMM;
        if(DEBUG_MEMORY)
        {
            Rprintf("\n### Object construction complete. Getting ready to rumble... ###\n\n");
        }
        if(DEBUG)
        {
            Rprintf("\n### Fitting HMM model parameters ###\n\n");
        }

        int verbose = INTEGER(sexpverbose)[0];

// memory allocation
        double** alpha = NULL;
        double** beta = NULL;
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
        allocateMemEmissionProb(&emissionProb, maxLen, K);

        int **flags = NULL;
        int *state2flag = NULL;
        RGETFLAGS(sexpflags, sexpstate2flag, &flags, &state2flag, nsample, T, K);

        int* couples = (int*)malloc(sizeof(int)*K);
        for(i=0; i<K; i++)
        {
            couples[i] = i;
        }

        double logLik = 0;
        for(n=0; n<nsample; n++)
        {
            if(fixedEmission == NULL)
            {
                myHMM->calcEmissionProbs(obs, emissionProb, T, n, flags, state2flag, NULL, isNaN, ncores, verbose, couples);
//Rprintf("%f\n", emissionProb[0][0]);
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
            for(t=0; t<T[n]; t++)
            {
                if(isNaN[n][t] == 0)
                {
                    logLik = logLik + log(c[t]);
                }
            }
        }

        deallocateMemAlpha(alpha, maxLen, K);
        deallocateMemBeta(beta, maxLen, K);
        deallocateMemRescFac(c, maxLen, K);
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
        free(start_d);
        free(couples);

        SEXP sexpLL;
        PROTECT(sexpLL = NEW_NUMERIC(1));
        NUMERIC_POINTER(sexpLL)[0] = -logLik;
        UNPROTECT(1);

        return sexpLL;

    }

}
