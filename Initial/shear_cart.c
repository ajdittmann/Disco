#include "../omega.h"
#include "../disco.h"

static double nu = 0.0;
static double gam = 0.0;
static double x0 = 0.0;
static double cs0 = 0.0;
static double sig0 = 0.0;
static double v0 = 0.0;
static double vx0 = 0.0;
static double vy0 = 0.0;
static double phi0_o_pi = 0.0;
static double rho_kappa = 0.0;
static int isothermal_flag = 0;
static int visc_flag = 0;
static double *t = NULL;

void get_xyz(double *, double *);
void get_vec_from_xyz(double *, double *, double *);
void get_vec_contravariant(double *, double *, double *);

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    nu = theDomain->theParList.viscosity;
    isothermal_flag = theDomain->theParList.isothermal_flag;
    visc_flag = theDomain->theParList.visc_flag;
    x0 = theDomain->theParList.initPar1;
    cs0 = theDomain->theParList.initPar2;
    sig0 = theDomain->theParList.initPar3;
    v0 = theDomain->theParList.initPar4;
    vx0 = theDomain->theParList.initPar5;
    vy0 = theDomain->theParList.initPar6;
    phi0_o_pi = theDomain->theParList.initPar7;
    rho_kappa = theDomain->theParList.initPar8;
    t = &(theDomain->t);
}

void initial(double *prim, double *x)
{
    double xyz[3];
    get_xyz(x, xyz);

    double cosp = cos(phi0_o_pi*M_PI);
    double sinp = sin(phi0_o_pi*M_PI);
    
    double X = cosp*xyz[0] + sinp*xyz[1];

    double rho, vx, vy, P;

    rho = exp(-rho_kappa * X);

    if(isothermal_flag)
        P = get_cs2(x) * rho / gam;
    else
        P = cs0*cs0 * rho / gam;

    double t0 = sig0*sig0 / (2 * nu);
    double sig2 = 2 * nu * (*t+t0);
    if(!visc_flag)
        sig2 = sig0*sig0;

    double xc = x0 + (vx0 + rho_kappa*nu)*(*t);

    vx = vx0;
    vy = v0 * exp(-0.5*(X-xc)*(X-xc) / sig2) * sig0/sqrt(sig2) + vy0;

    double Vxyz[3] = {cosp*vx - sinp*vy, sinp*vx + cosp*vy, 0};
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
