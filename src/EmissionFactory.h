
#ifndef EMISSIONFACTORY_HEADER
#define EMISSIONFACTORY_HEADER

#include "EmissionFunction.h"
#include <cmath>
#include "DebugConstants.h"
#include <list>
using namespace std;

class EmissionFactory
{
    public:
// create single Emission function
        virtual EmissionFunction* createEmissionFunction(ParamContainerEmissions *emissionParams, int parallel) = 0;
// create set of Emission functions
        virtual EmissionFunction* createEmissionFunctionMixed(list<EmissionFunction*> efb, ParamContainerEmissions *emissionParams){};
        ~EmissionFactory()
        {
            if(DEBUG_MEMORY)
            {
                printf("delete->EmissionFactory;\n");
            }
        }
};

EmissionFactory* createEmissionFactory(int whichone);
#endif
