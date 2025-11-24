#include "../disco.h"
#include "../omega.h"

static double gam = 0.0;
static double M = 0.0;
static double cx = 0.0;
static double cy = 0.0;
static double sig = 0.0;
static double cs0 = 0.0;
static double vx0 = 0.0;
static double vy0 = 0.0;
static double rho_atm = 0.0;

static int isothermal_flag = 0;

void get_xyz(const double *x, double *xyz);
void get_vec_from_xyz(double *x, double *vxyz, double *v);
void get_vec_contravariant(double *x, double *v, double *vc);
    
void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    M = theDomain->theParList.initPar1;
    cx = theDomain->theParList.initPar2;
    cy = theDomain->theParList.initPar3;
    sig = theDomain->theParList.initPar4;
    vx0 = theDomain->theParList.initPar5;
    vy0 = theDomain->theParList.initPar6;
    rho_atm = theDomain->theParList.initPar7;
    cs0 = theDomain->theParList.initPar8;
    isothermal_flag = theDomain->theParList.isothermal_flag;
}

void initial( double * prim , double * x )
{
    double rho0 = M * 2 / (M_PI*sig*sig * (1.0 - 4/(M_PI*M_PI)));

    double xyz[3];
    get_xyz(x, xyz);
    
    double r = sqrt((xyz[0]-cx)*(xyz[0]-cx) + (xyz[1]-cy)*(xyz[1]-cy));
    double R = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);


    double rho, vx, vy, P;

    double falloff_pow = 8.0;

    if(r < sig)
        rho = 0.5*rho0*(1 + cos(M_PI*r/sig)) + rho_atm;
    else if (R < 1.0)
        rho = rho_atm;
    else 
        rho = rho_atm / sqrt(1 + pow(R-1, 2*falloff_pow));

    if(rho > rho_atm)
    {
        double fac = 1 - pow(1-rho/rho0, 4);
        //fac = 1.0;
        vx = vx0 * fac;
        vy = vy0 * fac;
    }
    else
    {
        vx = 0.0;
        vy = 0.0;
    }
    double cs = cs0 * pow(rho/rho0, 0.5*(gam-1));
    P = rho * cs*cs/gam;

    if(isothermal_flag)
    {
        double cs2 = get_cs2(x);
        P = rho * cs2 / gam;
    }


    double Vxyz[3] = {vx, vy, 0};
    double V[3];
    get_vec_from_xyz(x, Vxyz, V);
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
