#include "../paul.h"
#include "../omega.h"
#include "../geometry.h"

static double gam  = 0.0;
static double beta  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;
static double eps = 0.0;
static int isothermal_flag = 0;

//double get_nu( const double *, const double *);
//double get_cs2( const double *, const double *);


void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
   isothermal_flag = theDomain->theParList.isothermal_flag;
   eps = theDomain->theParList.grav_eps;
   beta = theDomain->theParList.coolPar1;
}

void initial( double * prim , double * x ){

   double r = x[0];

   double omega2 = 1.0/(r*r*r); // - 1./(gam*Mach*Mach))/(r*r*r);
   double omega = sqrt(omega2);

   double rho = 1.0;
   double Pp = rho*(get_cs2(x)/gam + (gam-1)*9*nu*beta*omega/4);

   double visc = get_nu(x, prim);

   double Vrpz[3] = {-1.5*visc/r, r*omega, 0.0};
   double V[3];
   get_vec_from_rpz(x, Vrpz, V);
   get_vec_contravariant(x, V, V);

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = V[0];
   prim[UPP] = V[1];
   prim[UZZ] = V[2];


}
