
#ifndef MULTINOMIALFACTORY_HEADER
#define MULTINOMIALFACTORY_HEADER

#include "EmissionFactory.h"
#include "Multinomial.h"
#include "ParamContainerEmissions.h"
#include "DebugConstants.h"

class MultinomialFactory : public EmissionFactory
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
                return new Multinomial(emissionParams);
            }
            return NULL;

        }
        MultinomialFactory() {}
        ~MultinomialFactory() { }
};
#endif
