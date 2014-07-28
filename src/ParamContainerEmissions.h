
#ifndef PARAMCONTAINEREMISSIONS_HEADER
#define PARAMCONTAINEREMISSIONS_HEADER

#include "matUtils.h"
#include "MemoryAllocation.h"
#include "DebugConstants.h"
#include "EmissionFunction.h"
#include <cmath>

#include <list>
using namespace std;

#define MULTIVARIATEGAUSSIAN 1
#define BERNOULLI 2
#define JOINTLYINDEPENDENT 3

class ParamContainerEmissions
{
protected:
    // general
    int whichone;
    int D;
    int nsample;
    int* start;
    int* T;
    
    // gauss
    double **mu;
    double **sigma;
    double **inverseSigma; 
    double determinant;
    double regularize;
    double** gammaAux;

    // Bernoulli
    double p; 
    // more to come...
    
public:

    // one constructor for each emission function.
    int getWhichOne();

    ParamContainerEmissions(double **mu, double **sigma, double regularize, int D, int* start);
    ParamContainerEmissions(double p, int D, int* start);
    ParamContainerEmissions(int D);
    double** getGaussianMU();
    double** getGaussianSIGMA();
    double** getGaussianINVSIGMA();
    double getGaussianDET();
    double getGaussianREG();
    double** getGammaAux();
    int* getT();

    void setGaussianMU(double **mu);
    void setGaussianMUelement(double val, int d);
    void setGaussianSIGMA(double **sigma);
    void setGaussianSIGMAelement(double val, int d1, int d2);
    void setGaussianINVSIGMAelement(double val, int d1, int d2);
    void setGaussianDET(double val);
    void setDataVars(int nsample, int* T);
    void setDataVars(double** wrapper_gamma);

    void setGammaAux(double val, int n, int t);
    int getNsample();

    double getBernoulliP();
    void setBernoulliPelement(double val, int d);
    void setBernoulliP(double p);

    int getD();
    int* getStart();
    virtual ~ParamContainerEmissions();

};
ParamContainerEmissions** allocateParamContainerVector(int d);
#endif
