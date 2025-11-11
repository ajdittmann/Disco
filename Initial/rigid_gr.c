#include "../paul.h"

static double gamma_law = 0.0;

void setICparams( struct domain * theDomain ){
    gamma_law = theDomain->theParList.Adiabatic_Index;
}

void initial( double * prim , double * x ){
   double r = x[0];
   double phi = x[1];

   double rho = 1.0;
   double omega = 1.0;
   double P0 = 1.0e-2;


   double xc = 0.0;
   double yc = 0.0;
   double dx = r*cos(phi) - xc;
   double dy = r*sin(phi) - yc;
   double dr2 = dx*dx + dy*dy;
   double dr = sqrt(dr2);

   double Pp  = (gamma_law-1.0)/gamma_law * ((rho+gamma_law/(gamma_law-1.0)*P0)
                     * pow(1.0-dr2*omega*omega,-0.5*gamma_law/(gamma_law-1.0)) 
                     - rho);
   double U = sqrt(dr2)*omega / sqrt(1-dr2*omega*omega);

   double Ux = -dy/dr * U;
   double Uy =  dx/dr * U;

   double UR =  cos(phi) * Ux + sin(phi) * Uy;
   double UP = -sin(phi) * Ux + cos(phi) * Uy;

   double X = 0.0; 
   if( dx > 0.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = UR;
   prim[UPP] = r * UP;
   prim[UZZ] = 0.0;

   if(NUM_C > 5)
   {
       int q;
       for(q=5; q<NUM_C; q++)
           prim[q] = 0.0;
   }

   if( NUM_N>0 ) prim[NUM_C] = X;
}
