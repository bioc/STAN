
#ifndef BERNOULLIFACTORY_HEADER
#define BERNOULLIFACTORY_HEADER

#include "EmissionFactory.h"
#include "Bernoulli.h"
#include "ParamContainerEmissions.h"
#include "DebugConstants.h"

class BernoulliFactory : public EmissionFactory
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
                return new Bernoulli(emissionParams);
            }
            else if(parallel == 1)
            {
                return new Bernoulli(emissionParams);
            }
        }
        BernoulliFactory() {}
        ~BernoulliFactory() { }
};
#endif
