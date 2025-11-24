#include "../disco.h"

static double gam  = 0.0;
static double visc   = 0.0;
static double Mach = 0.0;
static double beta = 0.0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   visc   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
}

void initial( double * prim , double * x ){

   double mdot = 1.0;
   double pi3 = 3.0*M_PI;
   double r = x[0];
   double omega02 = 1.0/pow(r,3.);
   double omegaP2 = 0.0;
   double beta = 1.0/Mach;
   
   double omega = sqrt( omega02 - omegaP2 );

   double nu = visc*beta*beta*omega*r*r;
   double rho = mdot/(pi3*nu);
   double Pp = rho*beta*beta*r*r*omega*omega;
   //double Pp = rho/(r*Mach*Mach);

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = -1.5*nu/r;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;

}
