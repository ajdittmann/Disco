#include "../paul.h"

// Cylindrical relativistic blast wave in a uniform magnetic field,
// cf. Del Zanna et al 2007
// Del Zanna values are:
//      r0 = 1.0
//      B0 = 0.1
//      relFac = 1.0

static double gamma_law = 0.0;
static double r0 = 0.0;
static double B0 = 0.0;
static double relfac = 0.0;
static double x0 = 0.0;

void setICparams( struct domain * theDomain ){
    gamma_law = theDomain->theParList.Adiabatic_Index;
    r0 = theDomain->theParList.initPar1;
    B0 = theDomain->theParList.initPar2;
    relfac = theDomain->theParList.initPar3;
    x0 = theDomain->theParList.initPar4;
}

void initial( double * prim , double * x ){

   double r   = x[0];
   double phi = x[1];
   double R;
   double y0 = 0.0;

   double xx = r*cos(phi);
   double yy = r*sin(phi);
   R = sqrt(xx*xx + (yy-x0)*(yy-x0));

   double rhoL = 1.0e-2;
   double PL = relfac*relfac*1.0;
   double rhoR = 1.0e-4;
   double PR = relfac*relfac*5.0e-4;

   double rho = R > r0 ? rhoR : rhoL;
   double Pp  = R > r0 ? PR : PL;


   prim[RHO] = rho; 
   prim[PPP] = Pp;
   prim[URR] = 0.0; 
   prim[UPP] = 0.0;
   prim[UZZ] = 0.0; 

   if(NUM_C == 8)
   {
       prim[BRR] =  relfac * B0 * cos(phi);
       prim[BPP] = -relfac * B0 * sin(phi);
       prim[BZZ] = 0.0; 
   }

   double X = cos(phi)>0 ? 1.0 : 0.0;
   if( NUM_N>0) prim[NUM_C] = X;

}
