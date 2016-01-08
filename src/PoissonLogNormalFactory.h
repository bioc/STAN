
#ifndef POISSONLOGNAORMALFACTORY_HEADER
#define POISSONLOGNAORMALFACTORY_HEADER

#include "EmissionFactory.h"
#include "PoissonLogNormal.h"
#include "ParamContainerEmissions.h"
#include "DebugConstants.h"

class PoissonLogNormalFactory : public EmissionFactory
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
                return new PoissonLogNormal(emissionParams);
            }
            else if(parallel == 1)
            {
                return new PoissonLogNormal(emissionParams);
            }
        }
        PoissonLogNormalFactory() {}
        ~PoissonLogNormalFactory() { }
};
#endif
