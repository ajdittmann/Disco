#include "../disco.h"
#include "../omega.h"
#include "../geometry.h"

static double gam  = 0.0;
static double nu  = 1.0;
static double Mach = 0.0;
static double eps = 0.0;
static double p = 0.0;
static double q = 0.0;
static int threeD = 1;
static double chi = 0.0;
static double d = 0.0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   Mach = theDomain->theParList.Disk_Mach;
   eps = theDomain->theParList.grav_eps;

   q = theDomain->theParList.Cs2_Par;  //q
   p = theDomain->theParList.visc_par; //p
   nu = theDomain->theParList.viscosity;

   chi = 0.5*(3 - q);
   d = chi + p;

   if(theDomain->Nz > 1)
   {
      threeD = 1;
   }

}

void initial( double * prim , double * x ){

   double rpz[3];
   get_rpz(x, rpz);
   double R = rpz[0];
   double z = rpz[2];

   double Nu = nu*pow(fmax(x[0],1e-10), p);
   double cs2 = get_cs2(x);

   double rho = 1.0*pow(R, -p); // actually Sigma
   double Vrpz[3];
   if(threeD)
   {
     double H = sqrt( cs2*CUBE(R) );
     rho /= H*sqrt(2*M_PI); // rho_mid
     rho *= exp(-0.5*SQR(z/H));
     double omega2 = 1/CUBE(R);
     omega2 *= 1.0 - (d + q)*SQR(H/R) - 0.5*q*SQR(z/R);
     double omega = sqrt(omega2);
     Vrpz[0] = -1.5*Nu/R;
     Vrpz[1] = R*omega;
     Vrpz[2] = 0.0;
   }
   else
   {
     double omega2 = (1.0/(R*R*R))*(1.0 - 1.0/(Mach*Mach*gam)*(p+q)*pow(R, 1-q)  );
     double omega = sqrt(omega2);
     Vrpz[0] = -1.5*Nu/R;
     Vrpz[1] = R*omega;
     Vrpz[2] = 0.0;
   }

   double Pp = rho*cs2/gam;
   double V[3];
   get_vec_from_rpz(x, Vrpz, V);
   get_vec_contravariant(x, V, V);

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = V[0];
   prim[UPP] = V[1];
   prim[UZZ] = V[2];
}
