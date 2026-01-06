#include "../disco.h"

static int n = 0;
static double amp = 0.0;
static double width = 0.0;
static double P0 = 0.0;
static double alpha = 0.0;

static double gam = 0.0;

void setICparams(struct domain *theDomain)
{
    n = theDomain->theParList.initPar0;
    amp = theDomain->theParList.initPar1;
    width = theDomain->theParList.initPar2;
    alpha = theDomain->theParList.initPar3;
    P0 = theDomain->theParList.initPar4;

    gam = theDomain->theParList.Adiabatic_Index;
}

void initial(double *prim, const double *x)
{
    double r = x[0];
    double phi = x[1];
    
    double R = 1.0;
    double Om = 1.0;
    double rho0 = 1.0;

    double rho = rho0 * pow(r/R, alpha);
    double P = P0 + rho0*Om*Om*R*R*pow(r/R, alpha+2) / (alpha+2);

    double X, omega;
    
    if(r < R)
    {
        X = 0.0;
        omega = -Om;
    }
    else
    {
        X = 1.0;
        omega = Om;
    }

    double vr = 0.0;

    if(fabs(r-R) < width)
    {
        double f = 0.5 * (cos((r-R)*M_PI/width) + 1.0);
        vr = amp * cos(n*phi) * f;
    }
    
    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = vr;
    prim[UPP] = omega;
    prim[UZZ] = 0.0;
    
    if(NUM_N > 0)
    {
        prim[NUM_C] = X;
    }
}
