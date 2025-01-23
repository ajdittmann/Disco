#include "../omega.h"
#include "../paul.h"

/*
 * Disk-like flow in hydrostatic equilibrium. Setup assumes barotropic gas
 * (P = P(rho)) eg. ideal gas w/ uniform entropy (polytropic) or isothermal
 * gas with uniform temperature.
 *
 * Flow is in von Zeipel cylinders: omega = omega(r_cyl).
 */

static double gam  = 0.0;
static int isothermal = 0;
static int Hchoice = 0;
static double H0   = 0.0;
static double rho_atm0   = 0.0;
static double inv_beta   = 0.0;
static double bumpSig = 0.0;
static double bumpA = 0.0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   isothermal = theDomain->theParList.isothermal_flag;
   Hchoice = theDomain->theParList.initPar0;
   H0 = theDomain->theParList.initPar1;
   rho_atm0 = theDomain->theParList.initPar2;
   inv_beta = theDomain->theParList.initPar3;
   bumpA = theDomain->theParList.initPar4;
   bumpSig = theDomain->theParList.initPar5;
}

void get_rpz(double *x, double *rpz);
void get_coords_from_rpz(double *rpz, double *x);
void get_vec_from_rpz(double *x, double *Vrpz, double *V);
void get_vec_contravariant(double *x, double *V, double *Vc);

void initial( double * prim , double * X )
{
    double rpz[3];
    get_rpz(X, rpz);
    double r = rpz[0];
    double phi = rpz[1];
    double z = rpz[2];
    
    double GM = 1.0;
    double r0 = 1.0;

    double H, dHdr;

    if(Hchoice == 1)
    {
        H = H0 * r/r0;
        dHdr = H0/r0;
    }
    else if(Hchoice == 2)
    {
        H = H0 * sqrt(r*r*r/(r0*r0*r0));
        dHdr = 1.5 * H/r;
    }
    else if(Hchoice == 3)
    {
        //constant midplane density
        double hcoGM = 1.0/r0 - 1.0/sqrt(r0*r0+H0*H0);
        H = sqrt(2*r*r*r*hcoGM * (1-0.5*r*hcoGM) / ((1-r*hcoGM)*(1-r*hcoGM)));
        dHdr = H * (1.5/r - 0.25*hcoGM/(1-0.5*r*hcoGM) + hcoGM/(1-r*hcoGM));
    }
    else if(Hchoice == 4)
    {
        double rin = 0.5 * r0;
        double al = 1.5;

        double rpow = pow(r/r0, 2*al-1);

        H = H0 * sqrt((r-rin)/(r0-rin) * (rpow + 1));
        dHdr = H0*H0/(H*(r0-rin)) * 0.5*(rpow + 1 + (r-rin)/r*(2*al-1)*rpow);

        if(r <= rin)
        {
            H = 0;
            dHdr = 0;
        }
    }
    else
    {
        H = H0;
        dHdr = 0.0;
    }


    double RH = sqrt(r*r + H*H);
    double dRHdr = (r + H*dHdr) / RH;

    double PhiH = -GM/RH;
    double Phi = -GM/sqrt(r*r + z*z);
    double dPhiHdr = GM/(RH*RH) * dRHdr; 

    double om = sqrt(dPhiHdr / r);

    //enthalpy
    double h = -Phi + PhiH;
    double h0 = GM/r0 - GM/sqrt(r0*r0 + H0*H0);

    double rho, rho0, cs2, cs20, rho_atm, cs2_atm;

    double fac = 1.0;

    if(bumpSig > 0 && bumpA > 0)
    {
        double x = r*cos(phi);
        double y = r*sin(phi);
        double phi0 = 0.25 * M_PI;
        double x0 = 1.0 * cos(phi0);
        double y0 = 1.0 * cos(phi0);
        double d = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));

        if(d < bumpSig)
            fac += bumpA * 0.5*(1+cos(M_PI*d/bumpSig));
    }
    h *= fac;

    if(isothermal)
    {
        double rpz0[3] = {r0, rpz[1], 0.0};
        double x0[3];
        get_coords_from_rpz(rpz0, x0);

        cs20 = get_cs2(x0, 1.0);
        cs2 = get_cs2(X, 1.0);
        rho = rho_atm0 * exp(h/cs2);
        rho0 = rho_atm0 * exp(h0/cs20);

        rho_atm = rho_atm0;
        cs2_atm = cs2;
    }
    else
    {
        // rho0, cs20 are density, sound speed at r=r0, z=0
        rho0 = 1.0;

        cs20 = (gam-1.0) * h0;
        cs2 = (gam-1.0) * h;
        rho = rho0 * pow(cs2/cs20, 1.0/(gam-1.0));

        /*
        double h_atm = -Phi;
        cs2_atm = (gam-1) * h_atm;
        rho_atm = rho_atm0 * pow(cs2_atm/cs20, 1.0/(gam-1.0));
        */

        double beta = 6.0;
        double R = sqrt(r*r + z*z);
        cs2_atm = gam/beta * GM/R;
        rho_atm = rho_atm0 * pow(R/r0, -beta+1);
    }

    if(fabs(z) > H || rho < rho_atm || rho*cs2 < rho_atm*cs2_atm || H <= 0)
    {
        rho = rho_atm;
        cs2 = cs2_atm;
        om = 0.0;
    }

    double Pp = rho * cs2 / gam;

    double Bz = sqrt(2.0*rho0*cs20 * inv_beta);
    double Brpz[3] = {0.0, 0.0, Bz};
    double B[3];
    get_vec_from_rpz(X, Brpz, B);

    double Vrpz[3] = {0.0, r*om, 0.0};
    double V[3];
    get_vec_from_rpz(X, Vrpz, V);
    get_vec_contravariant(X, V, V);

    double Xpass = 0.0; 
    if( cos(phi) > 0.0 ) Xpass = 1.0; 

    prim[RHO] = rho;
    prim[PPP] = Pp;
    prim[URR] = V[0];
    prim[UPP] = V[1];
    prim[UZZ] = V[2];
    if(NUM_C > 5)
    {
        prim[BRR] = B[0];
        prim[BPP] = B[1];
        prim[BZZ] = B[2];
    }
    if( NUM_N>0 ) prim[NUM_C] = Xpass;
}
