
#include "../paul.h"

static double q_planet = 1.0;
static double Mach = 1.0;
static double eps = 0.0;
static double e_planet = 0.0;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = 2; 
   q_planet = theDomain->theParList.Mass_Ratio;
   Mach = theDomain->theParList.Disk_Mach;
   eps = theDomain->theParList.grav_eps;
   e_planet = theDomain->theParList.Eccentricity;
}

int planet_motion_analytic( void ){
   return(0);
}

void initializePlanets( struct planet * thePlanets ){

   double a  = 1.0;
   double e  = e_planet;
   double R = a*(1.-e);
   double om = pow( a , -1.5 )*sqrt(1.-e*e)/(1.-e)/(1.-e);

   double q = q_planet;
   double mu = q/(1.+q);

   thePlanets[0].M     = 1.-mu; 
   thePlanets[0].vr    = 0.0; 
   thePlanets[0].omega = om; 
   thePlanets[0].vz    = 0.0; 
   thePlanets[0].r     = R*mu; 
   thePlanets[0].phi   = M_PI; 
   thePlanets[0].z     = 0; 
   thePlanets[0].eps   = eps;
   thePlanets[0].type  = PLSPLINE;

   thePlanets[1].M     = mu;
   thePlanets[1].vr    = 0.0; 
   thePlanets[1].omega = om;  
   thePlanets[1].vz    = 0.0; 
   thePlanets[1].r     = R*(1.-mu); 
   thePlanets[1].phi   = 0.0; 
   thePlanets[1].z     = 0.0; 
   thePlanets[1].eps   = eps;
   thePlanets[1].type  = PLSPLINE;
}

void movePlanets( struct planet * thePlanets , double t , double dt ){
   thePlanets[0].phi += thePlanets[0].omega*dt;
   thePlanets[1].phi += thePlanets[1].omega*dt;
}
