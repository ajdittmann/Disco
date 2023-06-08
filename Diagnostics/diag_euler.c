
#include "../paul.h"
#include "../planet.h"
#include "../geometry.h"

static double gamma_law = 0.0;

void setDiagParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
}

int num_diagnostics(void){
   return(8);
}

int num_snapshot_rz(void)
{
    return 0;
}

int num_snapshot_arr(void)
{
    return 0;
}

/* Generic Diagnostics for Euler*/

void get_diagnostics( double * x , double * prim , double * Qrz, 
                        struct domain * theDomain )
{
   double r = x[0];

   double rho = prim[RHO];
   double vr = prim[URR];
   double omega = prim[UPP];
   double vz = prim[UZZ];
   double Pp = prim[PPP];

   Qrz[0] = rho;
   Qrz[1] = 2.*M_PI*r*rho*vr;
   Qrz[2] = 2.*M_PI*r*rho*vz;
   double rFp = 0.0;
   if( theDomain->Npl > 1 ){
      struct planet * pl = theDomain->thePlanets+1;
      double xyz[3];
      get_xyz(x, xyz);
      double Fxyz[3];
      planetaryForce( pl , xyz, Fxyz);
      rFp = xyz[0] * Fxyz[1] - xyz[1] * Fxyz[0];
   }
   Qrz[3] = 2.*M_PI*r*rho*(rFp);
   Qrz[4] = 2.*M_PI*r*rho*( r*r*omega*vr );
   Qrz[5] = omega;
   Qrz[6] = 2*M_PI*r * rho*r*r*omega;
   Qrz[7] = Pp;
}


void get_inst_diagnostics( double * x , double * prim , double * Qrz, 
                        struct domain * theDomain )
{
    //Silence is golden.
}

void get_snapshot_rz(const double *x, const double *prim, double *Qrz, 
                        struct domain * theDomain )
{
    //Silence is golden.
}

void get_snapshot_arr(const double *x, const double *prim, double *Qarr, 
                        struct domain * theDomain )
{
    //Silence is golden.
}
