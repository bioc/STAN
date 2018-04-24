
#ifndef MULTGAUSSFACTORY_HEADER
#define MULTGAUSSFACTORY_HEADER

#include "EmissionFactory.h"
#include "MultivariateGaussian.h"
#include "ParamContainerEmissions.h"
#include "DebugConstants.h"

class MultivariateGaussianFactory : public EmissionFactory
{
    public:
        EmissionFunction* createEmissionFunction(ParamContainerEmissions *emissionParams, int parallel)
        {
            if(DEBUG_MEMORY)
            {
                Rprintf("factory->create():");
            }
            if(parallel == 0)
            {
                return new MultivariateGaussian(emissionParams);
            }
            return NULL;
        }
        MultivariateGaussianFactory() {}
        ~MultivariateGaussianFactory() { }
};
#endif
