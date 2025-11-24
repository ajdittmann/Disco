#include "../disco.h"
#include "../omega.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;
static double eps = 0.0;

double get_nu(double *, double *);

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
   eps = theDomain->theParList.grav_eps;
}

void initial( double * prim , double * x ){

   double r = x[0];
   double R;
   double u = r/eps;
   if (u < 0.5) R = eps/( 16.*u*u/3. - 48.*pow(u, 4.)*0.2 + 32.*pow(u, 5.)*0.2 - 14.*0.2);
   else if (u < 1.0) R = eps/(1./(15.*u) + 32.*u*u/3. - 16.*u*u*u + 48.*pow(u, 4.)*0.2 - 32.*pow(u, 5.)/15. - 16.*0.2);
   else R = r;

   double omega02 = 1.0/pow(R,3.);
   double omegaP2 = (1./(Mach*Mach*R*R*R));

   double omega2 = fmax( (omega02 - omegaP2), 0.0 );
   double omega = sqrt(omega2);

   double cs2 = 1.0/(R*Mach*Mach);
   double visc = get_nu(x, prim);
   double rho = 1.0;
   if (nu > 0.0) rho = rho/nu;

   rho = rho * exp( - pow(r/0.15, -30.0));

   double Pp = rho*cs2;
   double X = 0.0; 
   if( r*cos(x[1]) > 0.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = -1.5*visc/(R);
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
