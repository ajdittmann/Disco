#include "../paul.h"

// To replicate Mignone, Journal of Computational Physics 270 (2014), 784-814

static double gam = 0.0;
static double csa = 0.0;
static double csw = 0.0;

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    csa = theDomain->theParList.initPar1;
    csw = theDomain->theParList.initPar2;
}

void initial(double *prim, double *x)
{
    double r = x[0];
    double phi = x[1];
    double z = x[2];

    double R = sqrt(r*r + z*z);

    double rho, vR, P;

    if(R >= 1.0)
    {
        rho = 0.25;
        vR = 0.0;
        P = csa*csa*rho/gam;
    }
    else
    {
        vR = tanh(5*R);
        rho = 1.0 / (R*R*vR);
        P = pow(rho, gam) * csw*csw / gam;
    }

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = vR * r/R;
    prim[UPP] = 0.0;
    prim[UZZ] = vR * z/R;
}
