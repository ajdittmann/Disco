#include "../disco.h"
#include "../omega.h"

double *t = NULL;

void setICparams( struct domain * theDomain ){

    t = &(theDomain->t);
}

void initial( double * prim , double * x ){
   double r = x[0];
   double phi = x[1];

   double rho = 1.0;
   double omega = 1.0;

   //if( r2 < 0.01 ) Pp += 1.0;

   int m = 4;
   rho += 0.1 * sin(m*(phi-omega*(*t))) * exp(-(r-0.5)*(r-0.5)/(2*0.1*0.1));
   
   double Pp  = 0.01;


   double X = 0.0; 
   if( cos(phi) > 0.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = 0.0;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;
}
