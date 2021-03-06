
#include "../paul.h"

static double eps = 0.0;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = 1; 
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
   thePlanets[0].RK_dM = 0.0;
   thePlanets[0].dM = 0.0;
}

void movePlanets( struct planet * thePlanets , double t, double dt ){
   //Silence is golden.
}

void forcePlanets( struct planet * thePlanets , double dt ){
   //Silence is golden.
}

