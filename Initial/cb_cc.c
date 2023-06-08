
#include "../paul.h"
#include "../geometry.h"
#include "../omega.h"

static double gam = 0.0;
static double Mach = 0.0;
static double sig0 = 0.0;
static double d0 = 0.0;
static double Rcav = 0.0;
static double Rout = 0.0;
static double v0 = 0.0;
static int threeD = 0;
static double rho_atm = 0.0;

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    Mach = theDomain->theParList.Disk_Mach;

    sig0 = theDomain->theParList.initPar1;
    d0 = theDomain->theParList.initPar2;
    Rcav = theDomain->theParList.initPar3;
    Rout = theDomain->theParList.initPar4;
    v0 = theDomain->theParList.initPar5;
    rho_atm = theDomain->theParList.initPar6;

    if(theDomain->Nz > 1)
    {
        threeD = 1;
    }
}

void initial(double *prim, double *x)
{
    double rpz[3];
    get_rpz(x, rpz);
    double r = rpz[0];
    double phi = rpz[1];

    double cs2 = get_cs2(x);

    double f = 1.0;
    if(Rout > 0.0)
        f = 1 - 1/(1 + exp(-2*(r-Rout)));

    double Sig = sig0 * ((1-d0)*exp(-pow(Rcav/r, 12)) + d0) * f;

    double om0 = sqrt((1 - 1/(Mach*Mach)) / (r*r*r));
    double n = 4;
    double omega = pow(pow(om0, -n) + 1, -1/n);
    double vr = v0 * sin(phi) * r * exp(-pow(r/3.5, 6));

    double rho;

    if(threeD)
    {
        double z = rpz[2];
        double H = r / Mach;
        double rho0 = Sig / (sqrt(2*M_PI) * H);

        double Phi0 = -1.0 / r;
        double Phi = -1.0 / sqrt(r*r + z*z);

        rho = rho0 * exp(-(Phi - Phi0) / cs2);
        if(rho < rho_atm)
            rho = rho_atm;
        vr *= rho/rho0;
    }
    else
        rho = Sig;

    if(rho < rho_atm)
        rho = rho_atm;

    double P = rho*cs2/gam;


    double Vrpz[3] = {vr, r*omega, 0.0};
    double V[3];
    get_vec_from_rpz(x, Vrpz, V);
    get_vec_contravariant(x, V, V);

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = V[0];
    prim[UPP] = V[1];
    prim[UZZ] = V[2];

    int q;
    for(q = NUM_C; q < NUM_Q; q++)
        prim[q] = 0.0;
}
