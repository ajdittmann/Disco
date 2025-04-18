#include "../paul.h"
#include "../omega.h"

static double A  = 0.0;
static double r1 = 3.5;
static double r2 = 6.5;
static double r0 = 0.0;
static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;
static double eps = 0.0;

double get_nu( const double *, const double *);

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
   eps = theDomain->theParList.grav_eps;
   r0 =0.5*(r1+r2);
}

void initial( double * prim , double * x ){

  double s = x[0];
  double phi = x[1];
  double z = x[2];

  double r = sqrt(s*s + z*z);
  double omega = sqrt(1.0/(r*r*r));

  double lx = 0.0;
  if (r>r1){
    if (r<r2){
      lx = 0.5*A*(1.0+sin(M_PI*((r-r0)/(r2-r1)) ));
    }
  else{
    lx = A;
    }
  }
  double lz = sqrt(1-lx*lx);

  double z0 = r*lx*cos(phi);
  double s0 = r*sqrt(1-lx*lx*cos(phi)*cos(phi));

  double ds = s0-s;
  double dz = z0-z;
  double xi = sqrt(ds*ds+dz*dz);

  double visc = get_nu(x, prim);
  double sigma = 1.0/visc;
  double cs2 = get_cs2(x, sigma);
  double h = sqrt(2.0*cs2)/omega

  double rho = sigma*exp(-xi*xi/(h*h));
  double Pp = rho*cs2/gam;

  prim[RHO] = rho;
  prim[PPP] = Pp;
  prim[URR] = -1.5*visc/(r);



  prim[UPP] = omega*lz;
  prim[UZZ] = omega*lx*r;


   double X = 0.0;
   if( cos(x[1]) > 0.0 ) X = 1.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
