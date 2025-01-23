
#include "../paul.h"
#include "../geometry.h"
#include "../omega.h"

static double gam = 0.0;
static int isothermal_flag = 0;
static double rho_ref = 1.0;
static double cs_ref = 0.0;
static double v0_ref = 0.0;
static double L = 0.0;

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    isothermal_flag = theDomain->theParList.isothermal_flag;
    cs_ref = theDomain->theParList.initPar1;
    v0_ref = theDomain->theParList.initPar2;
    L = theDomain->theParList.initPar3;
}

double poly_step(double x, double a, double b, double ya, double yb)
{
    double y;
    if(x <= a)
        y = ya;
    else if(x >= b)
        y = yb;
    else
    {
        double X = (2*x - (a+b)) / (b-a);
        double Y = 15.0*X/8.0 * (1.0 - 2.0*X*X/3.0 + X*X*X*X/5.0);
        y = 0.5*(yb-ya) * Y + 0.5*(ya+yb);
    }

    return y;
}

void initial(double *prim, double *x)
{
    double xyz[3], rpz[3];
    get_xyz(x, xyz);
    get_rpz(x, rpz);

    double rho = rho_ref;
    double P;
    if(isothermal_flag)
        P = get_cs2(x, rho) * rho / gam;
    else
        P = cs_ref*cs_ref * rho / gam;

    double r = rpz[0];
    double v = poly_step(r, L, 2*L, 0.0, v0_ref);
    double Vxyz[3] = {v, 0.0, 0.0};

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
    if(NUM_Q > 5)
        prim[5] = xyz[1] > 0 ? 1.0 : 0.0;
}
