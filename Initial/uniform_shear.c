
#include "../paul.h"

static double gam = 0.0;
static double rho = 1.0;
static double cs = 1.0;
static double vx = 0.0;
static double vy = 0.0;
static double vz = 0.0;
static double Bx = 0.0;
static double By = 0.0;
static double Bz = 0.0;
static double *t = NULL;

void get_xyz(const double *x, double *xyz);
void get_vec_from_xyz(double *x, double *vxyz, double *v);
void get_vec_contravariant(double *x, double *v, double *vc);
    
void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    vx = theDomain->theParList.initPar1;
    vy = theDomain->theParList.initPar2;
    vz = theDomain->theParList.initPar3;
    Bx = theDomain->theParList.initPar4;
    By = theDomain->theParList.initPar5;
    Bz = theDomain->theParList.initPar6;
    t = &(theDomain->t);
}

void initial( double * prim , double * x )
{
    double P = rho * cs*cs/gam;

    double Vxyz[3] = {vx, vy, vz};
    double Bxyz[3] = {Bx, By, Bz};
    double V[3], B[3];
    get_vec_from_xyz(x, Vxyz, V);
    get_vec_contravariant(x, V, V);
    get_vec_from_xyz(x, Bxyz, B);

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = V[0];
    prim[UPP] = V[1];
    prim[UZZ] = V[2];

    if(NUM_C == 8)
    {
        prim[BRR] = B[0];
        prim[BPP] = B[1];
        prim[BZZ] = B[2];
    }

    int q;
    for(q = NUM_C; q < NUM_Q; q++)
        prim[q] = 0.0;

    if(NUM_N > 0)
    {
        double xyz[3];
        get_xyz(x, xyz);
        if(xyz[0] < Vxyz[0]*(*t)-0.1)
            prim[NUM_C] = 1.0;
    }
}
