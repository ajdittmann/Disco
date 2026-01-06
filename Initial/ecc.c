#include "../disco.h"

static int visc_flag = 0;
static double nu = 0.0;
static double r0 = 0.0;
static double dr = 0.0;
static double dr2 = 0.0;
static int prof = 0;
static double sig0 = 0.0;
static double mach = 0.0;
static double gam = 0.0;
static double M = 0.0;
static double a = 0.0;
static double q = 0.0;

static double ecc = 0.0;

void setICparams( struct domain * theDomain ){
    visc_flag = theDomain->theParList.visc_flag;
    nu = theDomain->theParList.viscosity;
    prof = theDomain->theParList.initPar0;
    r0 = theDomain->theParList.initPar1;
    dr = theDomain->theParList.initPar2;
    dr2 = theDomain->theParList.initPar3;
    sig0 = theDomain->theParList.initPar4;
    ecc = theDomain->theParList.initPar5;
    mach = theDomain->theParList.Disk_Mach;
    gam = theDomain->theParList.Adiabatic_Index;
    
    M = 1.0;
    q = theDomain->theParList.Mass_Ratio;
    if(strcmp(PLANETS, "cm") == 0)
        a = 1.0;
    else
        a = 0.0;
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

double dpoly_step(double x, double a, double b, double ya, double yb)
{
    double dy;
    if(x <= a)
        dy = 0.0;
    else if(x >= b)
        dy = 0.0;
    else
    {
        double X = (2*x - (a+b)) / (b-a);
        double dX = 2.0/(b-a);
        double dY = 15.0/8.0 * (1.0 - 2.0*X*X + X*X*X*X) * dX;
        dy = 0.5*(yb-ya) * dY;
    }

    return dy;
}

void initial( double * prim , double * x ){

   double r = x[0];
   double phi = x[1];
   double l = r*(1+ecc*cos(phi));

   double rho, drhodr;
   double lm = r0-dr;
   double lp = r0+dr;

   double sig1 = 1.0;

   
   double lm2 = lm - dr2;
   double lp2 = lp + dr2;
   if(l > lm && l < lp)
   {
       rho = sig1;
       drhodr = 0.0;
   }
   else if(l>lm2 && l < lm)
   {
       rho = poly_step(l, lm2, lm, sig0, sig1);
       drhodr = dpoly_step(l, lm2, lm, sig0, sig1);
   }
   else if(l>lp && l < lp2)
   {
       rho = poly_step(l, lp, lp2, sig1, sig0);
       drhodr = dpoly_step(l, lp, lp2, sig1, sig0);
   }
   else
   {
       rho = sig0;
       drhodr = 0.0;
   }

   double cs20 = 1.0/r0 / (mach*mach);
   double Pp0 = sig1 * cs20 / gam;

   double Pp, dPdr;
   if(prof > 0)
   {
      Pp = Pp0 * pow(rho / sig1, gam);
      dPdr = (gam*Pp/rho) * drhodr;
   }
   else
   {
       Pp = Pp0;
       dPdr = 0.0;
   }

   //Orbital velocities for semi-latus rectum 'l' and eccentricity 'ecc'.
   double omega = sqrt(M*(1+ecc*cos(phi))/(r*r*r));
   double vr = sqrt(M/l) * ecc * sin(phi);

   double X = 0.0; 
   if( l>lm && l<lp ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;
}
