
#ifndef RACCESS_HEADER
#define RACCESS_HEADER

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

SEXP getListElement(SEXP list, const char *str);
#endif
