/**
 * @file    RWrapper.h
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
 * This file contains the function definition of the wrapper function for interfacing C++ code from R.
 */

#ifndef RWRAPPER_HEADER
#define RWRAPPER_HEADER

#include <R.h>
#include <Rdefines.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Error.h>
#include "MemoryAllocation.h"
#include "ParamContainerEmissions.h"
#include <Rembedded.h>
#include <R_ext/Rdynload.h>

#ifndef SUPPORT_OPENMP
    #define RUNOPENMPVERSION 0
#else
    #define RUNOPENMPVERSION 1
#endif


/**
 * Wrapper function for interfacing C++ code from R
 *
 * @return returns results to R.
 */
SEXP RHMMFit(SEXP sexpobs, SEXP sexppi, SEXP sexpA, SEXP sexpemission, SEXP sexptype, SEXP sexpdim, SEXP sexpregularize, SEXP sexpk, SEXP sexpmaxIters, SEXP sexpparallel, SEXP sexpflags, SEXP sexpstate2flag, SEXP sexpcouples, SEXP sexprevop, SEXP sexpverbose, SEXP sexpupdateTransMat, SEXP sexpfixedEmission, SEXP bidirOptimParams, SEXP emissionPrior, SEXP sexpeffectivezero, SEXP sepconvergence, SEXP sexpincrementalEM);
SEXP RHMMVITERBI(SEXP sexpobs, SEXP sexppi, SEXP sexpA,  SEXP sexpemission, SEXP sexptype, SEXP sexpdim, SEXP sexpk, SEXP sexpverbose, SEXP sexpfixedEmission);
SEXP RGETPOSTERIOR(SEXP sexpobs, SEXP sexppi, SEXP sexpA, SEXP sexpemission, SEXP sexptype, SEXP sexpdim, SEXP sexpk, SEXP sexpverbose, SEXP sexpfixedEmission, SEXP sexpncores, SEXP sexpflags, SEXP sexpstate2flag);
#endif

