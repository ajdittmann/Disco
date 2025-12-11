#include "../disco.h"
#include "../Calc/calc.h"

static double Npl = 3;
static double q_planet = 0.0;
static double q_bin = 1.0;
static double e_planet = 0.0;
static double e_bin = 0.0;
static double f_planet = 0.0;
static double a_planet = 3.0;
static double a_bin = 1.0;
static double eps = 0.05;
static double eps_planet = 0.05;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = Npl;
   q_bin = theDomain->theParList.Mass_Ratio;
   e_bin = theDomain->theParList.Eccentricity;
   q_planet = theDomain->theParList.planetPar1;
   e_planet = theDomain->theParList.planetPar2;
   a_planet = theDomain->theParList.planetPar3;
   f_planet = theDomain->theParList.planetPar4*M_PI;
   eps_planet = theDomain->theParList.planetPar5;
   eps = theDomain->theParList.grav_eps;

}

double root0( double E , double e ){
   return( E  - e*sin(E) );
}

double root1( double E , double e ){
   return( 1. - e*cos(E) );
}

int planet_motion_analytic( void ){
   return(1);
}

void initializePlanets( struct planet * thePlanets ){

   double TOL = 1e-14;

   double M = f_planet;
   double e = e_planet;

   //Newton-Rapheson to solve M = E - e*sin(E)
   //Necessary because of arbitrary phase offset
   double E = M;  //Guess value for E is M.
   double ff = root0( E , e ) - M;
   while( fabs(ff) > TOL ){
      double dfdE = root1( E , e );
      double dE = -ff/dfdE;
      E += dE;
      ff = root0( E , e ) - M;
   }
   double f = e*a_planet;
   double b = sqrt( fabs(SQR(a_planet) - f*f) );
   double en = -0.5/a_planet;
   double l = b*sqrt(2*fabs(en));

   double x = a_planet*cos(E)-f;
   double y = b*sin(E);
   double R   = sqrt(x*x+y*y);
   double phi = atan2(y,x);

   double vr = sqrt( fabs( 2.*en + 2./R - l*l/R/R ) );
   if( y<0.0 ) vr *= -1.;
   double mu = q_planet/(1.+q_planet);

   thePlanets[2].M   = mu;
   thePlanets[2].r   = R*(1.-mu);
   thePlanets[2].phi = phi;
   thePlanets[2].omega = l/R/R;
   thePlanets[2].vr = vr*(1.-mu);

   thePlanets[2].eps   = eps_planet;
   thePlanets[2].type  = PLPOINTMASS;

   double M_bin = 1.0 - mu;
   double R_bin = R*mu;
   double omega_bin = l/R/R;
   double vr_bin = vr*mu;
   double phi_bin = phi + M_PI;

   double xc = R_bin*cos(phi_bin);
   double yc = R_bin*sin(phi_bin);
   double vxc = -omega_bin*R_bin*sin(phi_bin) + vr_bin*cos(phi_bin);
   double vyc =  omega_bin*R_bin*cos(phi_bin) + vr_bin*sin(phi_bin);

   double R_in = a_bin*(1.-e_bin);
   double om = pow( a_bin , -1.5 )*sqrt(1.-SQR(e_bin))/SQR(1.-e_bin);

   double q = q_bin;
   double mu_bin = q/(1.+q);

   double r0 = R_in*mu_bin;
   double r1 = R_in*(1-mu_bin);
   double phi0 = M_PI;
   double phi1 = 0.0;
   double omega0 = om;
   double omega1 = om;
   double vr0 = 0.0;
   double vr1 = 0.0;


   double x0 = r0*cos(phi0) + xc;
   double y0 = r0*sin(phi0) + yc;
   double vx0 = -om*r0*sin(phi0) + vr0*cos(phi0) + vxc;
   double vy0 =  om*r0*cos(phi0) + vr0*sin(phi0) + vyc;

   double x1 = r1*cos(phi1) + xc;
   double y1 = r1*sin(phi1) + yc;
   double vx1 = -om*r1*sin(phi1) + vr1*cos(phi1) + vxc;
   double vy1 =  om*r1*cos(phi1) + vr1*sin(phi1) + vyc;

   r0 = sqrt(x0*x0 + y0*y0);
   phi0 = atan2(y0, x0);
   vr0 = vx0*cos(phi0) + vy0*sin(phi0);
   omega0 = (-vx0*sin(phi0) + vy0*cos(phi0))/r0;

   r1 = sqrt(x1*x1 + y1*y1);
   phi1 = atan2(y1, x1);
   vr1    =   vx1*cos(phi1) + vy1*sin(phi1);
   omega1 = (-vx1*sin(phi1) + vy1*cos(phi1))/r1;

   thePlanets[0].M     = (1.0 - mu_bin)*M_bin;
   thePlanets[0].vr    = vr0;
   thePlanets[0].omega = omega0;
   thePlanets[0].vz    = 0.0;
   thePlanets[0].r     = r0;
   thePlanets[0].phi   = phi0;
   thePlanets[0].z     = 0.0;

   thePlanets[0].eps   = eps;
   thePlanets[0].type  = PLPOINTMASS;

   thePlanets[1].M     = mu_bin*M_bin;
   thePlanets[1].vr    = vr1;
   thePlanets[1].omega = omega1;
   thePlanets[1].vz    = 0.0;
   thePlanets[1].r     = r1;
   thePlanets[1].phi   = phi1;
   thePlanets[1].z     = 0.0;

   thePlanets[1].eps   = eps;
   thePlanets[1].type  = PLPOINTMASS;


}

void movePlanets( struct planet * thePlanets , double t , double dt ){

  movePlanetsRK( Npl, thePlanets, dt);

}
