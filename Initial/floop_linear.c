#include "../disco.h"

static int orientation = 0;
static double gam = 0.0;
static double rho_ref = 1.0;
static double cs_ref = 1.0;
static double Rring = 0.0;
static double DRring = 0.0;
static double X0 = 0.0;
static double Y0 = 0.0;
static double Z0 = 0.0;
static double vx = 0.0;
static double vy = 0.0;
static double vz = 0.0;
    
void get_xyz(double *, double *);
void get_vec_from_xyz(double *, double *, double *);
void get_vec_contravariant(double *, double *, double *);

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    orientation = theDomain->theParList.initPar0;
    Rring = theDomain->theParList.initPar1;
    DRring = theDomain->theParList.initPar2;
    X0 = theDomain->theParList.initPar3;
    Y0 = theDomain->theParList.initPar4;
    Z0 = theDomain->theParList.initPar5;
    vx = theDomain->theParList.initPar6;
    vy = theDomain->theParList.initPar7;
    vz = theDomain->theParList.initPar8;
}

void initial(double *prim, double *x)
{
    /*
     * Linear advection of field loop
     */

    double xyz[3];
    get_xyz(x, xyz);

    double B0 = 1.0e-10 * sqrt(rho_ref);

    double dx = xyz[0] - X0;
    double dy = xyz[1] - Y0;
    double dz = xyz[2] - Z0;

    double dr, dh;

    if(orientation == 1)
    {
        dr = sqrt(dy*dy+dz*dz) - Rring;
        dh = dx;
    }
    else if(orientation == 2)
    {
        dr = sqrt(dz*dz+dx*dx) - Rring;
        dh = dy;
    }
    else
    {
        dr = sqrt(dx*dx+dy*dy) - Rring;
        dh = dz;
    }

    double dR = sqrt(dr*dr + dh*dh);
    double Bp;
    if(dR < DRring)
        Bp = B0 * 0.5*(1.0+cos(M_PI*dR/DRring));
    else
        Bp = 0.0;

    double Bx, By, Bz;
    if(orientation == 1)
    {
        Bx = 0.0;
        By = -dz*Bp/(dr+Rring);
        Bz = dy*Bp/(dr+Rring);
    }
    else if(orientation == 2)
    {
        Bx = dz*Bp/(dr+Rring);
        By = 0.0;
        Bz = -dx*Bp/(dr+Rring);
    }
    else
    {
        Bx = -dy*Bp/(dr+Rring);
        By = dx*Bp/(dr+Rring);
        Bz = 0.0;
    }
    
    double rho = rho_ref;
    double P = rho_ref * cs_ref*cs_ref/gam - 0.5*Bp*Bp;
    double Vxyz[3] = {vx, vy, vz};
    double Bxyz[3] = {Bx, By, Bz};

    //V is in contravariant basis, B in orthonormal
    double V[3], B[3];
    get_vec_from_xyz(x, Vxyz, V);
    get_vec_contravariant(x, V, V);
    get_vec_from_xyz(x, Bxyz, B);

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = V[0];
    prim[UPP] = V[1];
    prim[UZZ] = V[2];

    if(NUM_C > BZZ)
    {
        prim[BRR] = B[0];
        prim[BPP] = B[1];
        prim[BZZ] = B[2];
    }

    int q;
    for(q = NUM_C; q < NUM_Q; q++)
    {
        if(q == NUM_C && dR < DRring)
            prim[q] = 1.0;
        prim[q] = 0.0;
    }
}
