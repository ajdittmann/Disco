
#include "../paul.h"

static double gamma_law = 0.0;

void setDiagParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
}

int num_diagnostics(void){
   //return(85);
   return(5);
}

int num_snapshot_rz(void)
{
    return 0;
}

int num_snapshot_arr(void)
{
    return 0;
}

/* Generic Diagnostics for 2D boundary layers. Only good for m<40 modes.*/

void get_diagnostics( double * x , double * prim , double * Qrz, 
                        struct domain * theDomain )
{
   double r = x[0];
   double phi = x[1];
   double z = x[2];

   double rho = prim[RHO];
   double vr = prim[URR];
   double omega = prim[UPP];
   double vz = prim[UZZ];
   double Pp = prim[PPP];
   double vp = r*omega;

   Qrz[0] = rho;
   Qrz[1] = rho*omega;
   Qrz[2] = rho*vp;
   Qrz[3] = rho*vr;
   Qrz[4] = rho*vp*vr;

   //int i;
   //double val = r*vr*sqrt(rho)/M_PI;
   //for (i=0; i<40; i++) {
   //   Qrz[i*2+4] = val*cos(phi*(i+1));
   //   Qrz[i*2+5] = val*sin(phi*(i+1));
   //}
}

void get_snapshot_rz(const double *x, const double *prim, double *Qrz, 
                        struct domain * theDomain )
{
    // Silence is golden.
}

void get_snapshot_arr(const double *x, const double *prim, double *Qarr, 
                        struct domain * theDomain )
{
    // Silence is golden.
}


