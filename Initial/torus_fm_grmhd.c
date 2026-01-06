#include "../disco.h"
#include "../Hydro/frame.h"

static double M = 0.0;
static double a = 0.0;
static double gam = 0.0;
static double Rin = 0.0;
static double ell = 0.0;
static double B0 = 0.0;
static double rho_atm = 0.0;
static double cs2_atm = 0.0;

void setICparams( struct domain * theDomain ){
   gam = theDomain->theParList.Adiabatic_Index;
   M = theDomain->theParList.metricPar2;
   a = theDomain->theParList.metricPar3;
   Rin = theDomain->theParList.initPar1;
   ell = theDomain->theParList.initPar2;
   B0 = theDomain->theParList.initPar3;
   rho_atm = theDomain->theParList.Density_Floor;
   cs2_atm = theDomain->theParList.Pressure_Floor;
}

void initial( double * prim , double * x ){

   double R   = x[0];
   double phi = x[1];
   double th  = x[2];
   //double R = sqrt(r*r + z*z);
   //double costh = z / R;
   //double sinth = r / R;
   double sinth = sin(th);
   double r = R*sinth;

   double u0d2 =  r*r*(R-2*M) / (r*r*R - ell*ell*(R-2*M));
   double u0d = -sqrt(r*r * (R-2*M) / (r*r*R - ell*ell*(R-2*M)));
   double u0din = -sqrt(Rin*Rin*(Rin-2*M) / (Rin*Rin*Rin - ell*ell*(Rin-2*M)));
   double upd = -ell * u0d;
   double urd = 0.0;
   double uthd = 0.0;

   double K = 1.0e-3;
   double rho = pow((gam-1.) * (fabs(u0din/u0d)-1.0) / (K*gam), 1./(gam-1.));
   double P = K * pow(rho, gam);

   // Transform from Schwarzschild coordinates to Spherical Kerr-Schild
   
   double l0 = u0d;
   double lR = -2*M/(R-2*M) * u0d + urd;
   double lp = upd;
   double lth = uthd;

   if(R < Rin || u0d2 <= 0.0 || fabs(u0din) < fabs(u0d))
   {
       double al, be[3], gam[9], U[4], bed[3];
       int i, j;
       frame_U(x, U);
       al = metric_lapse(x);
       metric_shift(x, be);
       metric_gam(x, gam);
       for(i=0; i<3;i++)
       {
           bed[i] = 0.0;
           for(j=0; j<3; j++)
               bed[i] += gam[3*i+j]*be[j];
       }
       lR = bed[0]*U[0] + gam[0]*U[1] + gam[1]*U[2] + gam[2]*U[3];
       lp = bed[1]*U[0] + gam[3]*U[1] + gam[4]*U[2] + gam[5]*U[3];
       lth = bed[2]*U[0] + gam[6]*U[1] + gam[7]*U[2] + gam[8]*U[3];
       rho = rho_atm;
       P = rho_atm * cs2_atm;
   }

   prim[RHO] = rho;
   prim[PPP] = P;
   prim[URR] = lR;
   prim[UPP] = lp;
   prim[UZZ] = lth;

   if(NUM_C > BZZ)
   {
       double Bz = B0;
        prim[BRR] = 0.0;
        prim[BPP] = 0.0;
        prim[BZZ] = Bz;
   }

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}

