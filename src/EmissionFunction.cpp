
#include "EmissionFunction.h"

void EmissionFunction::computeShared(EmissionFunction** myEmissions, int nStates) {}
void EmissionFunction::resetShared() {}

EmissionFunction** allocateEmissionFunctionVector(int d)
{
    EmissionFunction **vector = (EmissionFunction**)malloc(sizeof(EmissionFunction*)*d);
    if(vector == NULL)
    {
        error("Not enough memory!\n");
    }
    return vector;
}


ParamContainerEmissions* EmissionFunction::getParameter()
{
    return this->emissionParams;
}


std::list<EmissionFunction*>  EmissionFunction::getEmissionFunctionList()
{
    return this->ContainerList;
}


// EmissionFunction::~EmissionFunction(){
// 	printf("yo\n");
// 	//delete this->emissionParams;
// }
