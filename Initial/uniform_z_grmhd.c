#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

    double fac = 1.0e-3;
    double bfac = 0.1;
    double vfac = 0.0;
   prim[RHO] = 1.0;
   prim[PPP] = 1.0*fac*fac;
   prim[URR] = 0.0;
   prim[UPP] = 0.0;
   prim[UZZ] = fac*vfac;
   if( NUM_N>0 ) prim[NUM_C] = 0.0;

   if(NUM_C > 5)
   {
       prim[BRR] = 0.0;
       prim[BPP] = 0.0;
       prim[BZZ] = fac*bfac;
   }
}
