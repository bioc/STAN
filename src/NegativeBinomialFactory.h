
#ifndef NEGATIVEBINOMIALFACTORY_HEADER
#define NEGATIVEBINOMIALFACTORY_HEADER

#include "EmissionFactory.h"
#include "NegativeBinomial.h"
#include "ParamContainerEmissions.h"
#include "DebugConstants.h"

class NegativeBinomialFactory : public EmissionFactory
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
                return new NegativeBinomial(emissionParams);
            }
            else if(parallel == 1)
            {
                return new NegativeBinomial(emissionParams);
            }
        }
        NegativeBinomialFactory() {}
        ~NegativeBinomialFactory() { }
};
#endif
