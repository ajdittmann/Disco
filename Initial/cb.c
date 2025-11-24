#include "../disco.h"
#include "../omega.h"

static double gam = 0.0;
static double visc = 0.0;
static double Mdot = 0.0;
static double R0 = 0.0;
static double rho0 = 0.0;

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    visc = theDomain->theParList.viscosity;
    Mdot = theDomain->theParList.initPar1;
    R0 = theDomain->theParList.initPar2;
    rho0 = theDomain->theParList.initPar3;
}

void initial(double *prim, double *x)
{
    double r = x[0];
    double phi = x[1];

    double om;
    if(r>1.0)
        om = pow(r, -1.5);
    else
        om = 1.0;

    double nu = get_nu(x, prim);

    double rho = Mdot/(3*M_PI*nu);
    double cs2 = get_cs2(x);
    if(nu == 0.0)
        rho = Mdot/(3*M_PI * 0.1*cs2/om);
    double v = -1.5*nu/r;

    double fac = (atan(5*(r-R0)/R0) + 0.5*M_PI) / M_PI;

    rho = (1-fac)*rho0 + fac*rho;
    v = fac*v;

    double P = rho*cs2/gam;

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = v;
    prim[UPP] = om;
    prim[UZZ] = 0.0;

    int q;
    for(q = 5; q < NUM_Q; q++)
        prim[q] = 0.0;
}
