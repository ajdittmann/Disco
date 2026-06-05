#include "../disco.h"
#include "../omega.h"
#include "../geometry.h"

static double gam  = 0.0;
static double nu  = 1.0;
static double Mach = 0.0;
static double eps = 0.0;
static double mu = 0.0;
static double p = 0.0;
static double q = 0.0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   Mach = theDomain->theParList.Disk_Mach;
   eps = theDomain->theParList.grav_eps;

   q = theDomain->theParList.Cs2_Par;  //q
   p = theDomain->theParList.visc_par; //p
   nu = theDomain->theParList.viscosity;
   mu = theDomain->theParList.Mass_Ratio;
}

void initial( double * prim , double * x ){

   double R = x[0];

   double Nu = nu*pow(fmax(x[0],1e-10), p);

   //double depth = 1.0/(1.0 + 25.0*nu/(Mach*SQR(Mach*mu)));
   double depth = 1.0/(1.0 + 3*25.0*nu/(Mach*SQR(Mach*mu)));
   double width = pow( nu*SQR(SQR(SQR(Mach)))  , -0.2);

   double rho = pow(R, -p)*(1.0 - depth*exp(-SQR(SQR((R-1.0)/width))));
   double Pp = rho*get_cs2(x)/gam;


   double omega2 = (1.0/(R*R*R))*(1.0 - 1.0/(Mach*Mach*gam)*(p+q)*pow(R, 1-q)  );
   double omega = sqrt(omega2);

   double Vrpz[3] = {-1.5*Nu/R, R*omega, 0.0};
   double V[3];
   get_vec_from_rpz(x, Vrpz, V);
   get_vec_contravariant(x, V, V);

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = V[0];
   prim[UPP] = V[1];
   prim[UZZ] = V[2];
}
