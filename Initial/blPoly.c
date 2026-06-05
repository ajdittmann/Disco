#include "../disco.h"
#include "../geometry.h"
#include "../omega.h"

static double gam  = 0.0;
static double gm1  = 0.0;
static double pref  = 0.0;
static double Mach = 0.0;

static double rbl = 1.0;
static double dbl = 0.025;
static double om0 = 0.0;

static double r_inn = 0.0;
static double r_out = 0.0;
static double C_mid = 0.0;
static double C_inn = 0.0;
static double A = 0.0;
static double B = 0.0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   Mach = theDomain->theParList.Disk_Mach;
   //rbl = theDomain->theParList.initPar1;
   //dbl = theDomain->theParList.initPar2;
   //om0 = theDomain->theParList.initPar3;
   gm1 = gam - 1.0;
   pref = gam/(gm1*Mach*Mach);


   r_inn = rbl - dbl;
   r_out = rbl + dbl;
   double om_outer = 1.0/sqrt(CUBE(r_out));
   double om_delta = om_outer - om0;

   A = om0 + om_delta*(dbl - 1.0)/(2.0*dbl);
   B = om_delta*0.5/dbl;

   C_mid = 0.5*SQR(A*r_out) + (2./3)*A*B*CUBE(r_out) + 0.25*SQR(B*SQR(r_out)) + 1/r_out - pref;
   double smid = 0.5*SQR(A*r_inn) + (2./3)*A*B*CUBE(r_inn) + 0.25*SQR(B*SQR(r_inn)) + 1/r_inn - C_mid;
   C_inn = 0.5*SQR(om0*r_inn) + 1.0/r_inn - smid;

}

void initial( double * prim , double * x ){
   double rpz[3];
   get_rpz(x, rpz);

   double r = rpz[0];

   double omega = 1.0/sqrt(CUBE(r));
   double rho = 1.0;

   if (r <= r_out && r > r_inn) {
     omega = A + r*B;
     rho = pow((0.5*SQR(A*r) + (2./3)*A*B*CUBE(r) + 0.25*SQR(B*SQR(r)) + 1/r - C_mid)/pref, 1/gm1);
   }
   if (r <= r_inn) {
     omega = om0;
     rho = pow((0.5*SQR(om0*r) + 1/r - C_inn)/pref, 1/gm1);
   }

   double Pp = pow(rho, gam)/SQR(Mach);

   double X = 0.0;
   if( r > rbl - dbl ) X = (r-(rbl-dbl))/(2.0*dbl);
   if( r > rbl + dbl ) X = 1.0;

   double vr = 0.0;
   double Vrpz[3] = {vr, r*omega, 0.0};
   double V[3];
   get_vec_from_rpz(x, Vrpz, V);
   get_vec_contravariant(x, V, V);

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = V[0];
   prim[UPP] = V[1];
   prim[UZZ] = V[2];
   if( NUM_N>0 ) prim[NUM_C] = X;

}
