#include "../omega.h"
#include "../paul.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;
static double q_planet	= 0.0;
static double a		= 1.0;
static double rot_om	= 0.0;
static double mach_csd	= 0.0;
static double alpha_csd= 0.0;
static int profChoice = 0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
   q_planet	= theDomain->theParList.Mass_Ratio;
   rot_om	= theDomain->theParList.RotOmega;
   profChoice	= theDomain->theParList.initPar0;
   mach_csd	= theDomain->theParList.initPar1;
   alpha_csd	= theDomain->theParList.initPar2; 
}


void initial( double * prim , double * x ){

   double r		= x[0];
   double phi		= x[1];
   double rx		= r*cos(phi);
   double ry		= r*sin(phi);
   double R		= sqrt( (a-rx)*(a-rx) + ry*ry );

   double mu		= q_planet/(1.+q_planet);
   double omega		= sqrt(1./(R*R*R));

   double rho;
   double rho_0 = 1.0;

   if(profChoice == 1)
   {

       //double alpha_visc	= nu*Mach*Mach/(a*a*omega);
       //double K		= q_planet*q_planet*Mach*Mach*Mach*Mach*Mach/alpha_visc;
       double K		= pow(q_planet, 2)*pow(mach_csd, 5)/alpha_csd;

       double f_0		= 0.45;
       //double tau_sh	= 1.89 + 0.53/(q_planet*Mach*Mach*Mach);
       //double tau_r		= 0.3363585661*pow(fabs(1.5*Mach*((R/a) - 1.0)), 2.5);
       double tau_sh	= 1.89 + 0.53/(q_planet*pow(mach_csd, 3));
       double tau_r		= 0.3363585661*pow(fabs(1.5*mach_csd*((R/a) - 1.0)), 2.5);
       double f_r;
       if (tau_r < tau_sh){
        f_r	= f_0;
       }else{
        f_r	= f_0*sqrt(tau_sh / tau_r);
       }
       
       double nom		= f_r*K/(3.0*M_PI);
       double denom		= 1.0 + f_0*K/(3.0*M_PI);
       
       rho		= rho_0 * (1.0 - (nom/denom)*sqrt(a/R));
   }
   else
       rho = rho_0;

   if(q_planet <= 0.0)
       rho = rho_0;

   //double Pp		= rho/(gam*Mach*Mach);
   //double Pp		= rho/(gam*mach_csd*mach_csd);
   double  Pp = rho * get_cs2(x, rho) / gam;

   double X = 0.0; 
   if( r*cos(x[1]) > 0.0 ) X = 1.0; 
   
   //double theta		= atan2(ry, a-rx);
   //double xi		= M_PI - theta - phi;
   double dr_dPhi	= -a*sin(phi);
   double dphi_dPhi	= 1.0 - (a/r)*cos(phi);

   double vr_global	= dr_dPhi*(omega - sqrt(1./(a*a*a)));
   double vp_global	= dphi_dPhi*(omega - sqrt(1./(a*a*a))) + rot_om;
   double vp_local	= sqrt(mu/(r*r*r));
   
   //Fixing angular velocity profile 
   double r_in		= 0.005;
   double r_out		= 0.01;
   double vr_f, vp_f, slope_p, slope_r, slope_rho;
   double rho_in	= 1.1*rho;
   
   if (r < r_in)
   {
	vp_f	= vp_local;
	vr_f	= 0;
	rho	= rho_in;
   }
   else if (r_in < r && r < r_out)
   {
	slope_p	= (vp_global - vp_local)/(r_out - r_in);
	vp_f	= vp_local + slope_p*(r-r_in);

	slope_r	= (vr_global - 0)/(r_out - r_in);
	vr_f	= 0 + slope_r*(r-r_in);
	
	slope_rho = (rho - rho_in)/(r_out - r_in);
	rho	= rho_in + slope_rho*(r-r_in);
   }else{
	vp_f	= vp_global;
	vr_f	= vr_global;
	rho	= rho;
   }

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr_f; //-sin(xi) * (sqrt(1./(R)) - sqrt(1./(a)); //-1.5*nu/r;
   prim[UPP] = vp_f; //cos(xi) R * sqrt(1./(a*a*a));
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
