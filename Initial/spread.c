#include "../paul.h"

static double nu = 0.0;
static double l0 = 0.0;
static double amp0 = 0.0;
static double r0 = 0.0;
static double sig0 = 0.0;
static double Mach = 0.0;
static double gam = 0.0;

void setICparams( struct domain * theDomain ){
   nu = theDomain->theParList.viscosity;
   l0 = theDomain->theParList.initPar1;
   amp0 = theDomain->theParList.initPar2;
   r0 = theDomain->theParList.initPar3;
   sig0 = theDomain->theParList.initPar4;
   Mach = theDomain->theParList.Disk_Mach;
   gam = theDomain->theParList.Adiabatic_Index;
}

void initial( double * prim , double * x ){

   double r = x[0];
   double phi = x[1];

   double rho0 = 1.0;
   // double Mdot0 = 3*M_PI*nu*rho0;
   // double Jdot0 = Mdot0*j0;
   double j = sqrt(r);

   double f = exp(-0.5*(r-r0)*(r-r0)/(sig0*sig0));

   double rho = rho0 * (1 - l0/j) + amp0 * f;
   double cs2 = 1.0/(Mach*Mach);
   double Pp  = rho * cs2 / gam;

   double drhodr = 0.5*l0/(j*j*j) - (r-r0)/(sig0*sig0) * amp0 * f;

   double vr = -1.5*nu/r * (1 + 2*r*drhodr/rho);

   double X = 0.0; 
   if( cos(phi) > 0.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr;
   prim[UPP] = j / (r*r) + 0.5 * drhodr * cs2 * r / (gam * rho * j);
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;
}
