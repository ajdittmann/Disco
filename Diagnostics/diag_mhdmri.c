#include "../disco.h"

static double gamma_law = 0.0;

void setDiagParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
}

int num_diagnostics(void){
   return(18);
}

int num_inst_diagnostics(void){
   return(0);
}

void planetaryForce( struct planet * , double , double , double , double * , double * , double * , int );

/* MHD */

void get_diagnostics( double * x , double * prim , double * Qrz, 
                        struct domain * theDomain )
{
   double r = x[0];
   double phi = x[1];
   double z = x[2];

   double rho = prim[RHO];
   double Pp = prim[PPP];
   double vr = prim[URR];
   double vp = prim[UPP];
   double vz = prim[UZZ];
   double Br = prim[BRR];
   double Bp = prim[BPP];
   double Bz = prim[BZZ];

   double cs2 = gamma_law * Pp / rho;
   double B2 = Br*Br + Bp*Bp + Bz*Bz;

   double sr = rho*vr;
   double sp = r*r*rho*vp;
   double sz = rho*vz;
   double e = 0.5*rho*(vr*vr + r*r*vp*vp + vz*vz) + Pp/(gamma_law-1.0) + 0.5*B2;
   Qrz[0] = rho;
   Qrz[1] = Pp;
   Qrz[2] = vr;
   Qrz[3] = vp;
   Qrz[4] = vz;
   Qrz[5] = Br;
   Qrz[6] = Bp;
   Qrz[7] = Bz;
   Qrz[8] = rho;
   Qrz[9] = e;
   Qrz[10] = sr;
   Qrz[11] = sp;
   Qrz[12] = sz;
   Qrz[13] = rho*vr*vp;
   Qrz[14] = Br*Bp;
   Qrz[15] = B2;
   Qrz[16] = sqrt(cs2);
   Qrz[17] = sqrt(B2/rho);

}

void get_inst_diagnostics( double * x , double * prim , double * Qrz, 
                        struct domain * theDomain )
{
    //Silence is golden.
}
