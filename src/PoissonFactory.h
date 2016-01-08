
#ifndef POISSONFACTORY_HEADER
#define POISSONFACTORY_HEADER

#include "EmissionFactory.h"
#include "Poisson.h"
#include "ParamContainerEmissions.h"
#include "DebugConstants.h"

class PoissonFactory : public EmissionFactory
{
    public:
        EmissionFunction* createEmissionFunction(ParamContainerEmissions *emissionParams, int parallel)
        {
            if(DEBUG_MEMORY)
            {
                printf("factory->create():");
            }
            if(parallel == 0)
            {
                return new Poisson(emissionParams);
            }

        }
        PoissonFactory() {}
        ~PoissonFactory() { }
};
#endif
