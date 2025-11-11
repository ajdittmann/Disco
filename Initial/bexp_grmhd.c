#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){
   double r   = x[0];
   double phi = x[1];

   double rho0 = 1.0e-2;

   double rhoI = rho0;
   double PpI   = 1.0;
   double rhoO = 1.0e-2 * rho0;
   double PpO   = 3.0e-5;
   double r1 = 0.8;
   double r2 = 1.0;
   
   double b = 0.1;//1.0e-3; //0.1;
   double phi0 = 0.0; //M_PI/4.;

   double rho, Pp;

   if( r<r1 )
   {
       rho = rhoI;
       Pp  =  PpI;
   }
   else if(r  < r2)
   {
       rho = exp(((r2-r)*log(rhoI) + (r-r1)*log(rhoO)) / (r2-r1));
       Pp  = exp(((r2-r)*log( PpI) + (r-r1)*log( PpO)) / (r2-r1));
   }
   else
   {
       rho = rhoO;
       Pp  =  PpO;
   }

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = 0.0;
   prim[UPP] = 0.0;
   prim[UZZ] = 0.0;

   prim[BRR] =  b*cos(phi-phi0);
   prim[BPP] = -b*sin(phi-phi0);
   prim[BZZ] = 0.0;

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}
