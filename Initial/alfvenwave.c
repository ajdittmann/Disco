#include "../disco.h"

static int waveChoice = 0;
static double gam = 0.0;
static double rho_ref = 1.0;
static double cs_ref = 1.0;
static double v0 = 0.0;         // Fluid speed (in B0 direction)
static double cA = 0.0;         // Ambient Alfven speed =sqrt(B0^2/rho)
static double a = 0.0;          // Amplitude of B wave (in units of B0)
static double lam = 0.0;        // Wavelength of wave
static double costheta0 = 0.0;  // (co?)latitude of B0 direction. = B0z / B0
static double phi0_o_pi = 0.0;  // longitude of B0 direction.=atan2(B0y,B0x)/pi
static double phase0 = 0.0;     // phase of wave at infinity (in units of pi)
static double L = 0.0;          // width of wave packet
static double x0 = 0.0;         // Location of wave packet center
    
void get_xyz(double *, double *);
void get_vec_from_xyz(double *, double *, double *);
void get_vec_contravariant(double *, double *, double *);

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    waveChoice = theDomain->theParList.initPar0;
    if(waveChoice == 0)
    {
        v0 = theDomain->theParList.initPar1;
        cA = theDomain->theParList.initPar2;
        a = theDomain->theParList.initPar3;
        lam = theDomain->theParList.initPar4;
        phi0_o_pi = theDomain->theParList.initPar5;
        costheta0 = theDomain->theParList.initPar6;
        phase0 = theDomain->theParList.initPar7;  
    }
    else
    {
        v0 = theDomain->theParList.initPar1;
        cA = theDomain->theParList.initPar2;
        a = theDomain->theParList.initPar3;
        lam = theDomain->theParList.initPar4;
        phi0_o_pi = theDomain->theParList.initPar5;
        costheta0 = theDomain->theParList.initPar6;
        x0 = theDomain->theParList.initPar7;  
        L = theDomain->theParList.initPar8;
    }
}

void initial(double *prim, double *x)
{
    /*
     * Plane-parallel (circularly polarized?) Alfven wave
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

    double X = (xyz[0]*ix+xyz[1]*iy+xyz[2]*iz - x0);

    double rho = rho_ref;
    double P = rho_ref * cs_ref*cs_ref/gam;
    double Vxyz[3] = {v0*ix, v0*iy, v0*iz};

    double B0 = sqrt(rho_ref * cA*cA);
    double Bxyz[3] = {B0*ix, B0*iy, B0*iz};

    // Alfven wave phase
    double psi;
    if(waveChoice == 0)
        //Simple harmonic wave
        psi = M_PI*phase0 + 2*M_PI*X/lam;
    else
        // Coiled wave packet
        psi = 2*M_PI*L/lam * 0.5*(1.0+tanh(X/L));
    
    //Perpendicular (waving) component of magnetic field
    double cp = cos(psi);
    double sp = sin(psi);
    double Bp[3] = {a*B0*(cp*jx+sp*kx), a*B0*(cp*jy+sp*ky), a*B0*(cp*jz+sp*kz)};

    //Add waving components to B and V
    int q;
    for(q=0; q<3; q++)
    {
        Bxyz[q] += Bp[q];
        Vxyz[q] += -Bp[q]/sqrt(rho_ref);
    }

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

    for(q = NUM_C; q < NUM_Q; q++)
        prim[q] = 0.0;
}
