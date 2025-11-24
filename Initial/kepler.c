#include "../disco.h"
#include "../omega.h"

static double gam  = 0.0;
static double Mach = 0.0;
static int isothermal_flag = 0;

double get_nu( const double *, const double *);


void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   Mach = theDomain->theParList.Disk_Mach;
   isothermal_flag = theDomain->theParList.isothermal_flag;
}

void initial( double * prim , double * x ){

   double r = x[0];

   double omega02 = 1.0/pow(r, 3);
   double omegaP2 = 0.0;

   double omega2 = fmax( (omega02 - omegaP2), 0.0 );
   double omega = sqrt(omega2);


   double rho = 1.0;
   double Pp = rho/(gam * Mach * Mach);
   if(isothermal_flag)
       Pp = rho * get_cs2(x);

   double X = 0.0;
   if( r*cos(x[1]) > 0.0 ) X = 1.0;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = 0.0;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
