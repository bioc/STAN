
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
        if(parallel == 0)
        {
            return new JointlyIndependent(efb, emissionParams);
        }
        else if(parallel == 1)
        {
            return new JointlyIndependent(efb, emissionParams);

        }
    // create Emission Function Mixed; wird überschriebn
    }

    JointlyIndependentFactory() {}
    ~JointlyIndependentFactory() { }
};
#endif
