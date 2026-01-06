//http://www.aanda.org/articles/aa/pdf/2010/08/aa12443-09.pdf
#include "../disco.h"

static double gam = 0.0;
static int perturbType = 0.0;
static double mach = 0.0;
static double n = 0;

void setICparams( struct domain * theDomain ){
    gam = theDomain->theParList.Adiabatic_Index;
    perturbType = theDomain->theParList.initPar0;
    mach = theDomain->theParList.initPar1;
    n = theDomain->theParList.initPar2;
}

void initial( double * prim , double * x ){

   double r   = x[0];
   //double phi = x[1];
   double z   = x[2];

   double rho = 1.0;
   double omega = pow( r , -1.5 );
   double cs20 = 1.0/mach;  //Keplerian velocity at r=R0=1 is 1
   double Pp = cs20*rho/gam;

   double Bz = 0.05513/n;
   if(n <= 0.0)
       Bz = 0.0;

   if( r < 2. || r > 3. ) Bz = 0.0;

   double vr = 0.0;
   double vz = 0.0;

   if(perturbType == 1)
   {
       double v0 = 5e-4;
       vr = v0*sin(2*M_PI*n*z);
   }
       

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr;
   prim[UPP] = omega;
   prim[UZZ] = vz;

   prim[BRR] = 0.0;
   prim[BPP] = 0.0;
   prim[BZZ] = Bz;

   if( NUM_N>0 ) prim[NUM_C] = (r<2.0||r>3.0) ? 0.0 : 1.0;

}

