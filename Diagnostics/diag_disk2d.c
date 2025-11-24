
#include "../disco.h"

static double gamma_law = 0.0;

void setDiagParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
}

int num_diagnostics(void){
   return(5);
}

int num_inst_diagnostics(void){
   return(0);
}

void planetaryForce( struct planet * , double , double , double , double * , double * , double * , int );

/* 2D Disk */

void get_diagnostics( double * x , double * prim , double * Qrz, 
                        struct domain * theDomain )
{
   double r = x[0];
   double phi = x[1];

   double rho = prim[RHO];
   double vr = prim[URR];
   double omega = prim[UPP];

   Qrz[0] = rho;
   Qrz[1] = 2.*M_PI*r*rho*vr;
   double Fr,Fp,Fz;
   Fp = 0.0;
   if( theDomain->Npl > 1 ){
      struct planet * pl = theDomain->thePlanets+1;
      planetaryForce( pl , r , phi , 0.0 , &Fr , &Fp , &Fz , 0 );
   }
   Qrz[2] = 2.*M_PI*r*rho*(r*Fp);
   Qrz[3] = 2.*M_PI*r*rho*( r*r*omega*vr );
   Qrz[4] = omega;
}


void get_inst_diagnostics( double * x , double * prim , double * Qrz, 
                        struct domain * theDomain )
{
    //Silence is golden.
}
