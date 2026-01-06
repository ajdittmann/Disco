#include "../disco.h"
#include "../geometry.h"
#include "../planet.h"

static double gamma_law = 0.0;

void setDiagParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
}

int num_diagnostics(void){
   return(5);
}

int num_snapshot_rz(void)
{
   return(3*2*30);
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

   double rho = prim[RHO];
   double vr = prim[URR];
   double omega = prim[UPP];
   double Pp = prim[PPP];
   double vp = r*omega;

   Qrz[0] = rho;
   Qrz[1] = Pp;
   Qrz[2] = rho*vp;
   Qrz[3] = rho*vr;
   Qrz[4] = rho*vp*vr;
}

void get_snapshot_rz(const double *x, const double *prim, double *Qrz, 
                        struct domain * theDomain )
{
   double rpz[3];
   get_rpz(x, rpz);
   double phi = rpz[1];
   double rho = prim[RHO];
   double V[3] = {prim[URR], prim[UPP], prim[UZZ]};
   get_vec_covariant(x, V, V);
   double Vrpz[3];
   get_vec_rpz(x, V, Vrpz);
   double vr = Vrpz[0];
   double vp = Vrpz[1];

   double vals[3] = {rho, vr, vp};

   double cospn, sinpn;
   double cosp = cos(phi);
   double sinp = sin(phi);
   double cospmn1 = cosp; // cos((n-1)*x)
   double sinpmn1 = sinp; // sin((n-1)*x)
   double cospmn2 = 1.0; // cos((n-2)*x)
   double sinpmn2 = 0.0; // sin((n-2)*x)

   int i,n;
   for(i=0; i<3; i++){
      Qrz[2*i]   = cospmn2*vals[i];
      Qrz[2*i+1] = sinpmn2*vals[i];
   }
   for(i=0; i<3; i++){
      Qrz[6+2*i]   = cospmn1*vals[i];
      Qrz[6+2*i+1] = sinpmn1*vals[i];
   }
   for(n=2; n<30; n++){
      cospn = 2*cospmn1*cosp-cospmn2;
      sinpn = 2*sinpmn1*cosp-sinpmn2;
      for(i=0; i<3; i++){
         Qrz[6*n+2*i] = cospn*vals[i];
         Qrz[6*n+2*i+1] = sinpn*vals[i];
      }
      cospmn2 = cospmn1;
      sinpmn2 = sinpmn1;
      cospmn1 = cospn;
      sinpmn1 = sinpn;
   }
}

void get_snapshot_arr(const double *x, const double *prim, double *Qarr, 
                        struct domain * theDomain )
{
    // Silence is golden.
}


