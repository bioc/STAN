
#ifndef JOINTLYINDEPENDENTFACTORY_HEADER
#define JOINTLYINDEPENDENTFACTORY_HEADER

#include "EmissionFactory.h"
#include "EmissionFunction.h"                     // new added
#include "MultivariateGaussian.h"
#include "JointlyIndependent.h"
#include "ParamContainerEmissions.h"
#include "DebugConstants.h"
#include <list>
using namespace std;

class JointlyIndependentFactory : public EmissionFactory
{
    public:
        EmissionFunction* createEmissionFunction(ParamContainerEmissions *emissionParams, int parallel)
        {
            if(DEBUG_MEMORY)
            {
                printf("createEmissionFunction factory->create():");
            }
            return NULL;
        }

        EmissionFunction* createEmissionFunctionMixed(list<EmissionFunction*> efb, ParamContainerEmissions *emissionParams)
        {
            int parallel=0;
            if(DEBUG_MEMORY)
            {
                printf("createEmissionFunctionMixed factory->create():");
            }
            return new JointlyIndependent(efb, emissionParams);
        }

        JointlyIndependentFactory() {}
        ~JointlyIndependentFactory() { }
};
#endif
