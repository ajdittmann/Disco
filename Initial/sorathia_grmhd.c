//http://www.aanda.org/articles/aa/pdf/2010/08/aa12443-09.pdf
#include "../disco.h"

static double Mach = 0.0;
static double gam = 0.0;
static double M = 0.0;
static double r0 = 0.0;
static double r1 = 0.0;
static double r2 = 0.0;
static double H0 = 0.0;
static int useB = 0;

void setICparams( struct domain * theDomain ){
   Mach = theDomain->theParList.Disk_Mach;
   gam = theDomain->theParList.Adiabatic_Index;
   M = theDomain->theParList.metricPar2;
   useB = theDomain->theParList.initPar0; //Turn on B field.
   r0 = theDomain->theParList.initPar1; // Inner edge
   r1 = theDomain->theParList.initPar2; // Fiducial radius
   r2 = theDomain->theParList.initPar3; // Outer Edge
   H0 = theDomain->theParList.initPar4; // Scale Height
}

double I(double r){
   double v=0.0;
   if( r>1.2*r0 && r<3.8/4.0 * r2 ) v = 1.0;
   return(v);
}

double floor( double x ){
   return( (double)(int)x );
}

double Rs( double x ){
   return( 1.2*r0 + (floor((x-r0)/H0)-.5)*H0 );
}

void initial( double * prim , double * x ){

   double rnd1 = (double)rand()/(double)(RAND_MAX);
   double rnd2 = (double)rand()/(double)(RAND_MAX);

   double r   = x[0];
   double phi = x[1];
   double z   = x[2];

   double rho0 = 100.0;
   double rho = rho0;
   double omega = sqrt(M/(r*r*r));
   double cs2;  // actually P/rho*h
   if(Mach > 0.0)
      cs2 = 1.0/(Mach*Mach);
   else
      cs2 = M/(2*r1*r1*r1) * H0*H0;

   double Pp = (gam-1.0) * cs2 / (gam*(1-cs2)-1) * rho;

   double dvp = (2e-2*rnd1 - 1e-2) * omega;
   double vz = (2e-2*rnd2 - 1e-2) * sqrt(M/r0);
   double vr = 0.0;
   double vp = omega +dvp;

   if(I(r) < 1.0)
   {
       vp = omega;
       vr = 0.0;
       vz = 0.0;
   }

   double A0 = 0.37 * sqrt(M/r0) * sqrt(rho0 / 100.0);
   if(!useB)
        A0 = 0.0;
   double Bz = A0*sin(2.*M_PI*(r-r0)/H0)*I(r)*(r0/r) * sqrt(r0/Rs(r));

   //Equatorial metric
   double g00 = -1+2*M/r;
   double g0r = 2*M/r;
   double grr = 1+2*M/r;
   double gzz = 1;
   double gpp = r*r;
   double u0 = 1.0 / sqrt(-g00 - 2*vr*g0r - grr*vr*vr - gzz*vz*vz - gpp*vp*vp);
   double lr = u0*(g0r+grr*vr); 
   double lp = u0*gpp*vp; 
   double lz = u0*(gzz*vz); 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = lr;
   prim[UPP] = lp;
   prim[UZZ] = lz;

   prim[BRR] = 0.0;
   prim[BPP] = 0.0;
   prim[BZZ] = Bz;

   int q;
   for(q=NUM_C; q<NUM_N; q++)
       prim[q] = 0.0;

}

