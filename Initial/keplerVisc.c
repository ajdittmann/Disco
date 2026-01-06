#include "../disco.h"
#include "../omega.h"
#include "../geometry.h"

static double gam  = 0.0;
static double nu  = 1.0;
static double beta  = 0.0;
static double Mach = 0.0;
static double eps = 0.0;
static int isothermal_flag = 0;
static double viscPar = 0.0;
static int viscChoice = 0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   Mach = theDomain->theParList.Disk_Mach;
   isothermal_flag = theDomain->theParList.isothermal_flag;
   eps = theDomain->theParList.grav_eps;
   beta = theDomain->theParList.coolPar1;

   viscChoice = theDomain->theParList.visc_profile;
   viscPar = theDomain->theParList.visc_par;
   nu = theDomain->theParList.viscosity;
}

void initial( double * prim , double * x ){

   double rs = sqrt(x[0]*x[0] + eps*eps);
   double r = x[0];

   //double nu = get_nu(x, prim);
   double rho = 1.0; ///nu;
   double Pp = rho*get_cs2(x)/gam;

   double dlogrho = 0.0;
   if (viscChoice == 1) dlogrho = -0.5;
   if (viscChoice == 2) dlogrho = -1*viscPar;

   double omega2 = (1.0/(rs*rs*r))*(1.0 - 1.0/(Mach*Mach) + dlogrho/(Mach*Mach));
   double omega = sqrt(omega2);

   double Vrpz[3] = {-1.5*nu/rs, r*omega, 0.0};
   double V[3];
   get_vec_from_rpz(x, Vrpz, V);
   get_vec_contravariant(x, V, V);

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = V[0];
   prim[UPP] = V[1];
   prim[UZZ] = V[2];
}
