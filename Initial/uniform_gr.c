#include "../disco.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

    double fac = 1.0e-1;
    double bfac = 0.0;
    double vfac = 0.0;
   prim[RHO] = 1.0;
   prim[PPP] = 1.0*fac*fac;
   prim[URR] = fac * vfac * cos(x[1]);
   prim[UPP] = fac * vfac * -sin(x[1])/x[0];
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = 0.0;

   if(NUM_C > 5)
   {
       prim[BRR] =  cos(x[1]) * fac * bfac;
       prim[BPP] = -sin(x[1]) * fac * bfac;
       prim[BZZ] = 0.0;
   }
}
