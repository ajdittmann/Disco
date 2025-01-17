#include "../paul.h"
#include "../omega.h"

static double gam = 0.0;
static double visc = 0.0;
static double Mdot = 0.0;
static double R0 = 0.0;
static double xi = 0.0;
static double delta = 0.0;
static double rho0 = 0.0;
static struct planet *thePlanets = NULL;
static double Mach = 0.0;
static int Npl = 0;

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    visc = theDomain->theParList.viscosity;
    Mdot = theDomain->theParList.initPar1;
    R0 = theDomain->theParList.initPar2;
    delta = theDomain->theParList.initPar5;
    xi = theDomain->theParList.initPar5;
    rho0 = theDomain->theParList.initPar3;
    Mach = theDomain->theParList.Disk_Mach;
    thePlanets = theDomain->thePlanets;
    Npl = theDomain->Npl;
}

void initial(double *prim, double *x)
{
    double r = x[0];
    double R = r/R0;
    double phi = x[1];


    double rho, fact;

    double nu = get_nu(x, prim);
    rho = 1.0;
    double om = 1.0;
    if (r > 2.0){
      rho = 1.00/nu;
      om = pow(r,-1.5);
    }

    double cs2 = get_cs2(x, rho);
    double v = -1.5*nu/(r+0.0001);
    double P = rho*cs2/gam;
 
    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = v;
    prim[UPP] = om;
    prim[UZZ] = 0.0;

    int q2;
    for(q2 = 5; q2 < NUM_Q; q2++)
        prim[q2] = 0.0;
}
