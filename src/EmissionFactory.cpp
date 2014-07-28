#include "EmissionFactory.h"
#include "MultivariateGaussianFactory.h"
#include "BernoulliFactory.h"
#include "JointlyIndependentFactory.h"
#include "ParamContainerEmissions.h"
#include <list>
using namespace std;

EmissionFactory* createEmissionFactory(int whichone)
{
    if (DEBUG_MEMORY)
    {
        printf("new->EmissionFactory;\n");
    }
    if (MULTIVARIATEGAUSSIAN == whichone)
    {
        return new MultivariateGaussianFactory();
    }
    if (BERNOULLI == whichone)
    {
        return new BernoulliFactory();
    }
    else
    {
        return new JointlyIndependentFactory();
    }
}
