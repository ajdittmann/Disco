#include "../disco.h"
#include "../Calc/calc.h"

static double q_planet = 1.0;
static double Npl = 3;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = Npl;
   q_planet = theDomain->theParList.Mass_Ratio;

}

int planet_motion_analytic( void ){
   return(1);
}

void initializePlanets( struct planet * thePlanets ){

   double a  = 1.0;

   double q = 1.0; //q_planet;
   double mu = q/(1.+q);

   double om = pow( a , -1.5 );

   thePlanets[0].M     = (1-mu);
   thePlanets[0].vr    = 0.0;
   thePlanets[0].omega = om;
   thePlanets[0].vz    = 0.0;
   thePlanets[0].r     = a*mu;
   thePlanets[0].phi   = M_PI;
   thePlanets[0].z     = 0.0;
   thePlanets[0].eps   = 0.05; //eps
   thePlanets[0].type  = PLPOINTMASS;

   thePlanets[1].M     = mu;
   thePlanets[1].vr    = 0.0;
   thePlanets[1].omega = om;
   thePlanets[1].vz    = 0.0;
   thePlanets[1].r     = a*(1.-mu);
   thePlanets[1].phi   = 0.0;
   thePlanets[1].z     = 0.0;
   thePlanets[1].eps   = 0.05; //eps;
   thePlanets[1].type  = PLPOINTMASS;


   thePlanets[2].M     = q_planet;
   thePlanets[2].vr    = 0.0;
   thePlanets[2].omega = pow(2.5, -1.5);
   thePlanets[2].vz    = 0.0;
   thePlanets[2].r     = 2.5;
   thePlanets[2].phi   = 0.0;
   thePlanets[2].z     = 0.0;
   thePlanets[2].eps   = 0.1;
   thePlanets[2].type  = PLPOINTMASS;
}

void movePlanets( struct planet * thePlanets , double t , double dt ){

  movePlanetsRK( Npl, thePlanets, dt);

}
