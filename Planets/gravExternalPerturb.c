
#include "../paul.h"

static double eps = 0.0;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = 3; 
   eps = theDomain->theParList.grav_eps;

}

int planet_motion_analytic( void ){
   return(1);
}

void initializePlanets( struct planet * thePlanets ){
   thePlanets[0].M     = 1.0; 
   thePlanets[0].vr    = 0.0; 
   thePlanets[0].omega = 1.0; 
   thePlanets[0].r     = 0.0; 
   thePlanets[0].phi   = 0.0; 
   thePlanets[0].eps   = eps;
   thePlanets[0].type  = PLPOINTMASS;

   thePlanets[1].M     = 1.0e-2; 
   thePlanets[1].vr    = 0.0; 
   thePlanets[1].omega = 1.0; 
   thePlanets[1].r     = 0.0; 
   thePlanets[1].phi   = 0.5*M_PI; 
   thePlanets[1].eps   = 0.0;
   thePlanets[1].type  = PLUNIFORM;

   thePlanets[2].M     = 1.0e-2; 
   thePlanets[2].vr    = 0.0; 
   thePlanets[2].omega = 1.0; 
   thePlanets[2].r     = 0.0; 
   thePlanets[2].phi   = 0.0; 
   thePlanets[2].eps   = 0.0;
   thePlanets[2].type  = PLLINEARX;
}

void movePlanets( struct planet * thePlanets , double t, double dt ){
   //Silence is golden.
}

