//http://www.aanda.org/articles/aa/pdf/2010/08/aa12443-09.pdf
#include "../disco.h"

static double M = 0.0;
static double gam = 0.0;

void setICparams( struct domain * theDomain ){
   gam = theDomain->theParList.Adiabatic_Index;
   M = theDomain->theParList.metricPar2;
}

void initial( double * prim , double * x ){

   double rnd1 = (double)rand()/(double)(RAND_MAX);
   double rnd2 = (double)rand()/(double)(RAND_MAX);

   double r   = x[0];
   double phi = x[1];
   double z   = 0.0; //x[2];

   double r0 = 1.0e6;
   double r1 = 2*r0;
   double r2 = 3*r0;


   double n = 4.0;

   double rho = 1.0;
   double vp = sqrt(M / (r*r*r));
   double cs = 0.1*vp*r;
   double Pp = cs*cs*rho/gam;

   double vr = (1e-3*rnd1 - 5e-4) * sqrt(M/r0);
   double vz = (1e-3*rnd2 - 5e-4) * sqrt(M/r0);
   double Bz = 0.05513/n * sqrt(M/r0);

   if (r < 0.5*(3*r1-r2) || r > 0.5*(3*r2-r1))
   {
       vr = 0.0;
       vz = 0.0;
   }
   if (r < r1 || r > r2)
       Bz = 0.0;

   double R2 = r*r+z*z;
   double R = sqrt(R2);
   double R3 = R*R2;

   double g00 = -1+2*M/R;
   double g0r = 2*M*r/R2;
   double g0z = 2*M*z/R2;
   double grr = 1+2*M*r*r/R3;
   double grz = 2*M*r*z/R3;
   double gzz = 1+2*M*z*z/R3;
   double gpp = r*r;
   double u0 = 1.0 / sqrt(-g00 - 2*vr*g0r - 2*vz*g0z - grr*vr*vr - 2*grz*vr*vz
                            - gzz*vz*vz - gpp*vp*vp);
   double lr = u0*(g0r+grr*vr+grz*vz); 
   double lp = u0*gpp*vp; 
   double lz = u0*(g0z+grz*vr+gzz*vz); 


   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = lr;
   prim[UPP] = lp;
   prim[UZZ] = lz;

   prim[BRR] = 0.0;
   prim[BPP] = 0.0;
   prim[BZZ] = Bz;

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}

