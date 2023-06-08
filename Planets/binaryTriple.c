
#include "../paul.h"

static double q_binary = 1.0;
static double q_bary = 1.0;
static double a_binary = 1.0;
static double a_bary = 1.0;
static double Mach = 1.0;
static double eps = 0.0;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = 3;
   q_bary = theDomain->theParList.Mass_Ratio;
   q_binary = theDomain->theParList.initPar1;
   //a_bary = theDomain->theParList.RotD;
   a_binary = theDomain->theParList.initPar2;

   Mach = theDomain->theParList.Disk_Mach;
   eps = theDomain->theParList.grav_eps;
}

int planet_motion_analytic( void ){
   return(1);
}

void initializePlanets( struct planet * thePlanets ){

   double q = q_bary;
   double mu = q/(1.+q);

   double mu2 = q_binary/(1.+q_binary);

   double om = sqrt(mu/pow(a_binary,3.0));

   thePlanets[0].M     = 1.-mu;
   thePlanets[0].vr    = 0.0; 
   thePlanets[0].omega = 0.0; 
   thePlanets[0].vz    = 0.0; 
   thePlanets[0].r     = a_bary;
   thePlanets[0].phi   = 0.0;
   thePlanets[0].z     = 0.0;
   thePlanets[0].eps   = 0.0;
   thePlanets[0].type  = PLPOINTMASS;


   thePlanets[1].M     = mu*(1.-mu2);
   thePlanets[1].vr    = 0.0; 
   thePlanets[1].omega = om; 
   thePlanets[1].vz    = 0.0; 
   thePlanets[1].r     = a_binary*mu2; 
   thePlanets[1].phi   = M_PI; 
   thePlanets[1].z     = 0.0; 
   thePlanets[1].eps   = eps;
   thePlanets[1].type  = PLPOINTMASS;


   thePlanets[2].M     = mu*mu2;
   thePlanets[2].vr    = 0.0;
   thePlanets[2].omega = om;
   thePlanets[2].vz    = 0.0;
   thePlanets[2].r     = a_binary*(1.-mu2);
   thePlanets[2].phi   = 0.0;
   thePlanets[2].z     = 0.0;
   thePlanets[2].eps   = eps;
   thePlanets[2].type  = PLPOINTMASS;

}

void movePlanets( struct planet * thePlanets , double t , double dt ){
   thePlanets[1].phi += thePlanets[1].omega*dt;
   thePlanets[2].phi += thePlanets[2].omega*dt;
}

