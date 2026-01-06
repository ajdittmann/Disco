#include "../disco.h"
#include "../Calc/calc.h"

static double M = 0.0;
static double gam = 0.0;
static double rs = 0.0;
static double R0 = 0.0;
static double B0 = 0.0;

void get_rpz(double *x, double *rpz);
void get_vec_from_rpz(double *x, double *vrpz, double *v);
void get_vec_contravariant(double *x, double *v, double *vc);

void setICparams( struct domain * theDomain ){
   gam = theDomain->theParList.Adiabatic_Index;
   M = 1.0;
   rs = theDomain->theParList.initPar1;
   R0 = theDomain->theParList.initPar2;
   B0 = theDomain->theParList.initPar3;
}

void initial( double * prim , double * x ){

    double rpz[3];
    get_rpz(x, rpz);
    double r = rpz[0];
    double z = rpz[2];

    double R = sqrt(r*r+z*z);
    double Mdot = 1.0;
    
    double us2 = M / (2*rs);
    double as2 = us2;
    double us = -sqrt(us2);

    double rhos = -Mdot / (4*M_PI*rs*rs*us);
    double K = as2 / (gam * pow(rhos,gam-1));

    double a02 = 0.5*(5.0-3*gam)*as2;
    double rho0 = pow(a02 / (gam * K), 1.0/(gam-1));


    double uR = bondi_newt_solve(Mdot, M, R, gam, rho0, sqrt(a02), K);
    double rhoB = -Mdot / (4*M_PI*R*R*uR);
    double PB = K * pow(rhoB, gam);

    double chi = R>R0 ? 1.0 : 0.5*(1+tanh(tan(M_PI*(R/R0-0.5))));

    double rho = rhoB;
    double P = PB;

    if(R < R0)
    {
        rho = rhos;
        P = K * pow(rhos, gam);
        uR = -sqrt(2*M/R);
    }

    double vrpz[3] = {r/R*uR*chi, 0, z/R*uR*chi};
    double v[3];
    get_vec_from_rpz(x, vrpz, v);
    get_vec_contravariant(x, v, v);

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = v[0];
    prim[UPP] = v[1];
    prim[UZZ] = v[2];

    if(NUM_C > 5)
    {
        double bb = rs*rs*B0/(R*R);
        double Brpz[3] = {bb*r/R, 0.0, bb*z/R};
        double B[3];
        get_vec_from_rpz(x, Brpz, B);
        prim[BRR] = B[0];
        prim[BPP] = B[1];
        prim[BZZ] = B[2];
    }

    int q;
    for(q=NUM_C; q < NUM_Q; q++)
    {
        prim[q] = 0.0;
    }
}
