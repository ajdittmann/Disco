
#include "../disco.h"

static double q_planet = 1.0;
static double eps = 0.0;
static double cosI = 1.0;
static double sinI = 0.0;
static double mu;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = 2;
   q_planet = theDomain->theParList.Mass_Ratio;
   eps = theDomain->theParList.grav_eps;
   sinI = sin(M_PI*theDomain->theParList.planetPar1/180);
   cosI = cos(M_PI*theDomain->theParList.planetPar1/180);
   mu = q_planet/(1.0 + q_planet);

}

int planet_motion_analytic( void ){
   return(1);
}

void initializePlanets( struct planet * thePlanets ){
   double rp = mu;
   double px = rp;
   double py = 0.0;
   double pz = 0.0;

   thePlanets[0].r = sqrt(SQR(px)+SQR(py));
   thePlanets[0].phi = atan2(py, px);
   thePlanets[0].z = pz;
   thePlanets[0].vr = 0.0;
   thePlanets[0].omega = cosI;
   thePlanets[0].vz = rp*sinI;
   thePlanets[0].f = 0.0;

   thePlanets[0].M     = 1.-mu;
   thePlanets[0].eps   = eps;
   thePlanets[0].type  = PLPOINTMASS;

   rp = 1.0-mu;
   px = -rp;
   py = 0.0;
   pz = 0.0;
   thePlanets[1].r = sqrt(SQR(px)+SQR(py));
   thePlanets[1].phi = atan2(py, px);
   thePlanets[1].z = pz;
   thePlanets[1].vr = 0.0;
   thePlanets[1].omega = cosI;
   thePlanets[1].vz = rp*sinI;
   thePlanets[1].f = M_PI;

   thePlanets[1].M     = mu;
   thePlanets[1].eps   = eps;
   thePlanets[1].type  = PLPOINTMASS;
}

void movePlanets( struct planet * thePlanets , double t , double dt ){
   UNUSED(t);

   thePlanets[0].f += dt;
   thePlanets[1].f += dt;

   double f = thePlanets[0].f;
   double cosF = cos(f); double sinF = sin(f);
   double rp = mu;
   double px = rp*cosF;
   double py = rp*sinF*cosI;
   double pz = rp*sinF*sinI;
   thePlanets[0].r = sqrt(SQR(px)+SQR(py));
   thePlanets[0].phi = atan2(py, px);
   thePlanets[0].z = pz;
   thePlanets[0].vr = rp*cosF*sinF*(SQR(cosI)-1);
   thePlanets[0].vz = rp*cosF*sinI;

   f = thePlanets[1].f;
   cosF = cos(f); sinF = sin(f);

   rp = 1.0-mu;
   px = rp*cosF;
   py = rp*sinF*cosI;
   pz = rp*sinF*sinI;
   thePlanets[1].r = sqrt(SQR(px)+SQR(py));
   thePlanets[1].phi = atan2(py, px);
   thePlanets[1].z = pz;
   thePlanets[1].vr = rp*cosF*sinF*(SQR(cosI)-1);
   thePlanets[1].vz = rp*cosF*sinI;
}
