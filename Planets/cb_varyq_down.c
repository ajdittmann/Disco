
#include "../paul.h"

static double q_planet = 1.0;
static double Mach = 1.0;
static double eps = 0.0;
static double nu = 1.0;
static double tau = 2.0;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = 2;
   q_planet = theDomain->theParList.Mass_Ratio;
   Mach = theDomain->theParList.Disk_Mach;
   eps = theDomain->theParList.grav_eps;
   nu = theDomain->theParList.viscosity;
}

int planet_motion_analytic( void ){
   return(1);
}

void initializePlanets( struct planet * thePlanets ){

   double a  = 1.0;
 
   double q = q_planet;
   double mu = q/(1.+q);

   double om = pow( a , -1.5 );

   thePlanets[0].M     = 1.-mu;
   thePlanets[0].vr    = 0.0;
   thePlanets[0].omega = om;
   thePlanets[0].vz    = 0.0;
   thePlanets[0].r     = a*mu;
   thePlanets[0].phi   = M_PI;
   thePlanets[0].z     = 0.0;
   thePlanets[0].eps   = eps;
   thePlanets[0].type  = PLPOINTMASS;

   thePlanets[1].M     = mu;
   thePlanets[1].vr    = 0.0;
   thePlanets[1].omega = om;
   thePlanets[1].vz    = 0.0;
   thePlanets[1].r     = a*(1.-mu);
   thePlanets[1].phi   = 0.0;
   thePlanets[1].z     = 0.0;
   thePlanets[1].eps   = eps;
   thePlanets[1].type  = PLPOINTMASS;
}

void movePlanets( struct planet * thePlanets , double t , double dt ){
   thePlanets[0].phi += thePlanets[0].omega*dt;
   thePlanets[1].phi += thePlanets[1].omega*dt;

   double T = tau*2.0*M_PI/nu;
   double q = 1.0/(1.+pow(t/T, 2.));
   double mu = q/(1.+q);

   thePlanets[0].M = 1.-mu;
   thePlanets[1].M = mu;

   thePlanets[0].r = mu;
   thePlanets[1].r = 1.-mu;

}
