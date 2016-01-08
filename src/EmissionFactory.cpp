#include "EmissionFactory.h"
#include "MultivariateGaussianFactory.h"
#include "BernoulliFactory.h"
#include "PoissonFactory.h"
#include "MultinomialFactory.h"
#include "NegativeBinomialFactory.h"
#include "PoissonLogNormalFactory.h"
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
    else if (BERNOULLI == whichone)
    {
        return new BernoulliFactory();
    }
    else if(POISSON == whichone)
    {
        return new PoissonFactory();
    }
    else if(MULTINOMIAL == whichone)
    {
        return new MultinomialFactory();
    }
    else if(JOINTLYINDEPENDENT == whichone)
    {
        return new JointlyIndependentFactory();
    }
    else if(NEGATIVEBINOMIAL == whichone)
    {
        return new NegativeBinomialFactory();
    }
    else if(POISSONLOGNORMAL == whichone)
    {
        return new PoissonLogNormalFactory();
    }
    else
    {
        error("Cannot create unknown emission factory!");
    }

}
