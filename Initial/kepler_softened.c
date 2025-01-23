
#include "../paul.h"
#include "../omega.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;
static double eps = 0.0;
static int isothermal_flag = 0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
   isothermal_flag = theDomain->theParList.isothermal_flag;
   eps = 0.1;
}

void initial( double * prim , double * x ){

   double r = x[0];

   double rho = 1.0;
   double Pp;
   if(isothermal_flag)
   {
       double cs2 = get_cs2(x, rho);
       Pp = rho*cs2/gam;
   }
   else
       Pp = rho/(Mach*Mach*gam);

   double omega = sqrt(1.0 / pow(r*r + eps*eps, 1.5));

   double X = 0.0; 
   if( r*cos(x[1]) > 0.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = -1.5*nu/r;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
