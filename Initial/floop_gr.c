#include "../paul.h"

static double gamma_law = 0.0;
static double x0 = 0.0;
static double Rl = 0.0;

void setICparams( struct domain * theDomain ){
    gamma_law = theDomain->theParList.Adiabatic_Index;
    x0 = theDomain->theParList.initPar1;
    Rl = theDomain->theParList.initPar2;
}

void initial( double * prim , double * x ){

   double r   = x[0];
   double phi = x[1];

   double B0 = 1.0e-4;
   double Om = 1.0;
   double P0 = 1.0e-2;
   double rho0 = 1.0;


   double n = gamma_law/(gamma_law-1.0);
   double rhoh = rho0 + n*P0;
   double Pp = (gamma_law-1.0)/gamma_law * (rhoh * pow(1-r*r*Om*Om,-0.5*n)
                    - rho0);

   double xl = r*cos(phi)-x0;
   double yl = r*sin(phi);

   double rl = sqrt(xl*xl+yl*yl);
   double xx = M_PI*rl/Rl;

   double Bp = B0*pow(sin(xx),2.)*sqrt(2.*rl/Rl);
   if( xx > M_PI) Bp = 0.0; 

   //double dP = B0*B0*( - (rl/Rl)*pow(sin(xx),4.) - (1./16./M_PI)*( 12.*xx - 8.*sin(2.*xx) + sin(4.*xx) ) ); 

   double Bx = -Bp*yl/rl;
   double By =  Bp*xl/rl;

   prim[RHO] = rho0; 
   prim[PPP] = Pp;
   prim[URR] = 0.0; 
   prim[UPP] = r*r*Om/sqrt(1.0-r*r*Om*Om);
   prim[UZZ] = 0.0; 

   prim[BRR] =  Bx*cos(phi) + By*sin(phi);
   prim[BPP] = -Bx*sin(phi) + By*cos(phi);
   prim[BZZ] = 0.0; 

   double X = cos(phi)>0 ? 1.0 : 0.0;
   if( NUM_N>0) prim[NUM_C] = X;

}
