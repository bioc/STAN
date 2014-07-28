
#include "matUtils.h"
using namespace std;
void inverse(double** A, int N)
{
    double *Atemp = (double*)malloc(sizeof(double)*N*N);
    int i,j;
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            Atemp[j+i*N] = A[i][j];
        }
    }

    int *IPIV = new int[N+1];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;

    F77_NAME(dgetrf)(&N,&N,Atemp,&N,IPIV,&INFO);
    if(INFO != 0)
    {
        error("Error in LU-Decomposition of covariance matrix.\n");
    }
    F77_NAME(dgetri)(&N,Atemp,&N,IPIV,WORK,&LWORK,&INFO);
    if(INFO != 0)
    {
        error("Error inverting covariance matrix.\n");
    }

    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            A[i][j] = Atemp[j+i*N];
        }
    }

    free(Atemp);
    delete IPIV;
    delete WORK;
}


void inverseR(double** mat, int N)
{
    int i,j;

    SEXP sexpmat;
    PROTECT(sexpmat = NEW_NUMERIC(N*N));
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            REAL(sexpmat)[j+N*i] = mat[i][j];
        }
    }

    // call solnp from R for optimization
    SEXP call = PROTECT( lang2( install( "c2invertCOV"), sexpmat ) ) ;
    SEXP res = PROTECT( eval( call, R_GlobalEnv ) ) ;

    // write results into matrix
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            Rprintf("%f ", REAL(sexpmat)[j+N*i] * mat[i][j]);
            mat[i][j] = REAL(sexpmat)[j+N*i];
        }
        Rprintf("\n");
    }

    UNPROTECT(1);
}


void **matrixMult(double **v1, int d11, int d12, double **v2, int d21, int d22, double **result)
{
    if(d12 != d21)
    {
        error("Wrong dimensions for matrix multiplication!\n");
    }
    int i,j,k;

    for(i=0; i<d11; i++)
    {
        for(j=0; j<d22; j++)
        {
            result[i][j] = 0;
            for(k=0; k<d12; k++)
            {
                result[i][j] = result[i][j] + v1[i][k]*v2[k][j];
            }
        }
    }
}


double matrixDet(double **m, int dim)
{
    int myNCol = dim;

    double  *myAP = new double[myNCol*(myNCol + 1)/2],
    *myW = new double[myNCol],
    *myZ = new double[myNCol*myNCol],
    *myWork = new double[myNCol * 3] ;
    int myInfo,
        myN = (int)(myNCol),
        myldz = (int)(myNCol) ;

    for (register int i = 0 ; i < myN ; i++)
        for (register int j = i ; j < myldz ; j++)
            myAP[i+(j+1)*j/2]  = m[i][j] ;

    F77_NAME(dspev)("V", "U", &myN, myAP, myW, myZ, &myldz, myWork, &myInfo) ;

    if (myInfo != 0)
        error("Non inversible matrix") ;
    double theDet;
    theDet = 1.0L ;
    for (register int i = 0 ; i < myNCol ; i++)
    {
        theDet *= myW[i] ;
    }

    delete myAP;
    delete myW;
    delete myZ;
    delete myWork;

    return theDet;
}
