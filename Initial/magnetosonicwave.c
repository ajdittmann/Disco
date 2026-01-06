#include "../disco.h"
#include "../Calc/calc.h"

static double gam = 0.0;
static double rho_ref = 1.0;
static double cs_ref = 1.0;
static double v0 = 0.0;         // Fluid speed (in B0 direction)
static double cA0 = 0.0;         // Ambient Alfven speed =sqrt(B0^2/rho)
static double a = 0.0;          // Amplitude of B wave (in units of B0)
static double L = 0.0;          // width of wave packet
static double costheta0 = 0.0;  // (co?)latitude of B0 direction. = B0z / B0
static double phi0_o_pi = 0.0;  // longitude of B0 direction.=atan2(B0y,B0x)/pi
static double x0 = 0.0;         // Location of wave packet center
static double phiB_o_pi = 0.0;     // Direction of B field w.r.t. wave
    
void get_xyz(double *, double *);
void get_vec_from_xyz(double *, double *, double *);
void get_vec_contravariant(double *, double *, double *);

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    v0 = theDomain->theParList.initPar1;
    cA0 = theDomain->theParList.initPar2;
    a = theDomain->theParList.initPar3;
    L = theDomain->theParList.initPar4;
    phi0_o_pi = theDomain->theParList.initPar5;
    costheta0 = theDomain->theParList.initPar6;
    x0 = theDomain->theParList.initPar7;
    phiB_o_pi = theDomain->theParList.initPar8;
}

void initial(double *prim, double *x)
{
    /*
     * Magnetosonic wave with purely transverse magnetic field
     */

    double xyz[3];
    get_xyz(x, xyz);

    double phi0 = phi0_o_pi * M_PI;
    double sintheta0 = sqrt((1.0-costheta0)*(1.0+costheta0));
    double cosphi0 = cos(phi0);
    double sinphi0 = sin(phi0);

    double ix = cosphi0*sintheta0;
    double iy = sinphi0*sintheta0;
    double iz = costheta0;

    double jx = -sinphi0;
    double jy = cosphi0;
    double jz = 0.0;

    double kx = -cosphi0*costheta0; //iy*jz - iz*jy
    double ky = -sinphi0*costheta0;//iz*jx - ix*jz
    double kz = sintheta0; //ix*jy - iy*jx

    double cpb = cos(phiB_o_pi * M_PI);
    double spb = sin(phiB_o_pi * M_PI);

    double bx = cpb*jx + spb*kx;
    double by = cpb*jy + spb*ky;
    double bz = cpb*jz + spb*kz;

    double X = (xyz[0]*ix+xyz[1]*iy+xyz[2]*iz - x0) / L;

    double rho = rho_ref;

    if(fabs(X) < 1.0)
    {
        double f = 1.0-X*X;
        rho = rho_ref * (1.0 + a*f*f*f*f);
    }

    double P = rho_ref*cs_ref*cs_ref/gam * pow(rho/rho_ref, gam);

    double v = v0 + magnetosonic_cf_int_newt(rho, rho_ref, cs_ref, 
                                                    cA0, gam);
    double Vxyz[3] = {v*ix, v*iy, v*iz};

    double B0 = sqrt(rho_ref)*cA0 * rho/rho_ref;
    double Bxyz[3] = {B0*bx, B0*by, B0*bz};

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
        prim[q] = 0.0;
}
