#include "../paul.h"
#include "../geometry.h"
#include "../omega.h"
#include "../planet.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;
static double eps = 0.0;
static int isothermal_flag = 0;

struct planet *thePlanets = NULL;
double R_cav = 0.0;
double R_out = 0.0;
double sig_atm = 0.0;


void setICparams( struct domain * theDomain )
{
    gam  = theDomain->theParList.Adiabatic_Index;
    nu   = theDomain->theParList.viscosity;
    Mach = theDomain->theParList.Disk_Mach;
    isothermal_flag = theDomain->theParList.isothermal_flag;
    eps = theDomain->theParList.grav_eps;
    thePlanets = theDomain->thePlanets;

    R_cav = theDomain->theParList.initPar1;
    R_out = theDomain->theParList.initPar2;
    sig_atm = theDomain->theParList.initPar3;
}


double psi(double x)
{
    if(x <= 1.4e-3)
        return 0.0;
    return exp(-1.0/x);
}

double psi_prime(double x)
{
    if(x <= 1.4e-3)
        return 0.0;
    return exp(-1.0/x) / (x*x);
}

double step0(double x)
{
    // step(x) = { 0 if x < 0, 1 if x > 1 and smoothly interpolates between}
    double a = psi(x);
    double b = psi(1-x);
    return a / (a + b);
}

double step0_prime(double x)
{
    // step(x) = { 0 if x < 0, 1 if x > 1 and smoothly interpolates between}
    double a = psi(x);
    double da = psi_prime(x);
    double b = psi(1-x);
    double db = -psi_prime(1-x);
    return (da*b - a*db) / ((a+b)*(a+b));
}

double smooth_interp(double x, double a, double b, double fx, double gx)
{
    double w = step0((x-a)/(b-a));
    return (1-w) * fx + w * gx;
}

void smooth_interp_d(double x, double a, double b, double fx, double gx,
                     double dfx, double dgx, double *f, double *df)
{
    double w = step0((x-a)/(b-a));
    double dw = step0_prime((x-a)/(b-a)) / (b-a);
    *f = (1-w) * fx + w * gx;
    *df = -dw * fx + (1-w) * dfx + dw * gx + w * dgx;
}

void initial( double * prim , double * x )
{
    double xyz[3];
    double rpz[3];
    get_xyz(x, xyz);
    get_rpz(x, rpz);

    //double logsig_inn = log(sig_atm);
    double logsig_inn = 0.5 * log(sig_atm);
    double logsig_mid = 0.0;
    double logsig_out = log(sig_atm);

    double logsig1 = 0.0;
    double dlogsig1 = 0.0;
    double logsig = 0.0;
    double dlogsig = 0.0;

    double r = rpz[0];

    //smooth_interp_d(r, 0.5*R_cav, R_cav, logsig_inn, logsig_mid, 0.0, 0.0,
    //                &logsig1, &dlogsig1);
    smooth_interp_d(r, 0.0, R_cav, logsig_inn, logsig_mid, 0.0, 0.0,
                    &logsig1, &dlogsig1);
    smooth_interp_d(r, R_out, 1.5*R_out, logsig1, logsig_out, dlogsig1, 0.0,
                    &logsig, &dlogsig);

    double sig = exp(logsig);
    double dsig = dlogsig * sig;

    double gxyz[3];
    planetaryForce(thePlanets, xyz, gxyz);

    double g[3], grpz[3];
    get_vec_from_xyz(x, gxyz, g);
    get_vec_rpz(x, g, grpz);

    double gr = grpz[0];

    double cs2 = get_cs2(x);
    double Pp = sig * cs2 / gam;

    double vr = 0.0;
    double om = sqrt((-gr*(1 - 1.0/(Mach*Mach)) + dsig*cs2/sig) / r);

    double vrpz[3] = {vr, r*om, 0.0};
    double v[3];
    get_vec_from_rpz(x, vrpz, v);
    get_vec_contravariant(x, v, v);

    prim[RHO] = sig;
    prim[PPP] = Pp;
    prim[URR] = v[0];
    prim[UPP] = v[1];
    prim[UZZ] = v[2];

    if(NUM_Q > NUM_C)
    {
        double X = (sig - sig_atm) / (1.0 - sig_atm);
        prim[NUM_C] = X;
    }
}
