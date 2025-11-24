#include "../disco.h"

static double gam = 0.0;
static double rho_ref = 1.0;
static double cs = 1.0;
static double mach = 0.0;
static double L = 0.0;
static double a = 0.0;
static double x0 = 0.0;
static double costheta0 = 0.0;
static double phi0_o_pi = 0.0;
static double *t = NULL;
    
void get_xyz(double *, double *);
void get_vec_from_xyz(double *, double *, double *);
void get_vec_contravariant(double *, double *, double *);

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    mach = theDomain->theParList.initPar1;
    a = theDomain->theParList.initPar2;     //1.0 in RAM
    x0 = theDomain->theParList.initPar3;  
    L = theDomain->theParList.initPar4;     //0.3 in RAM
    phi0_o_pi = theDomain->theParList.initPar5;
    costheta0 = theDomain->theParList.initPar6;

    t = &(theDomain->t);
}

void initial(double *prim, double *x)
{
    double xyz[3];
    get_xyz(x, xyz);

    double phi0 = phi0_o_pi * M_PI;
    double sintheta0 = sqrt(1.0-costheta0*costheta0);

    double kx = cos(phi0)*sintheta0;
    double ky = sin(phi0)*sintheta0;
    double kz = costheta0;

    double v0 = cs*mach;

    double X = (xyz[0]*kx+xyz[1]*ky+xyz[2]*kz - (x0 + v0*(*t))) / L;

    double rho, v, P;

    double f = X*X-1.0;
    if(fabs(X) < 1.0)
        rho = rho_ref*(1.0 + a*f*f*f*f);
    else
        rho = rho_ref;

    P = rho_ref * cs*cs/gam;

    v = v0;

    double Vxyz[3] = {v*kx, v*ky, v*kz};

    //Vxyz[0] = sin(4*xyz[0]) * cos(7*xyz[1]);
    //Vxyz[1] = cos(5*xyz[0] - 3*xyz[1]) - 1.0;

    double V[3];
    get_vec_from_xyz(x, Vxyz, V);
    get_vec_contravariant(x, V, V);

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = V[0];
    prim[UPP] = V[1];
    prim[UZZ] = V[2];

    int q;
    for(q = 5; q < NUM_Q; q++)
        prim[q] = 0.0;
}
