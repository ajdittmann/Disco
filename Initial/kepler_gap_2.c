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

struct Gap_Vals{
	double rho, d_r_rho;
};

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
   q_planet	= theDomain->theParList.Mass_Ratio;
   rot_om	= theDomain->theParList.RotOmega;
   mach_csd	= theDomain->theParList.initPar1;
   alpha_csd	= theDomain->theParList.initPar2; 
}

struct Gap_Vals gap_density(double r, double d, double M, double alpha, double q, double rho_0);
//double gap_density(double r, double d, double M, double alpha, double q, double rho_0);
//double d_gap_density(double r, double d, double M, double alpha, double q, double rho_0);
/*
double gap_density(double r, double a, double Mach, double alpha, double q, double rho_0){

	double qNL	= 1.04*pow(Mach, -3);
	double qw	= 34*qNL*pow(alpha*Mach, 0.5);
	double D	= 7*pow(Mach, 1.5)*pow(alpha, 0.25);
	double rho, A, del_q, denom, qr, q_denom;
	
	q_denom	= 1 + D*D*D*pow( pow( r/a, 1./6.) - 1, 6);
	qr	= q*pow(q_denom, -1./3.);
	
	if (qr < qNL){
		del_q = 1;
	} else { 
		del_q	= pow(qr/qNL, -0.5) + pow(qr/qw, 3);
	}
	
	A	= (0.45/(3*M_PI))*qr*qr*pow(Mach, 5)/alpha;
	denom	= 1 + A*del_q;
	
	return rho_0/denom;
}
*/
struct Gap_Vals gap_density(double r, double d, double M, double alpha, double q, double rho_0){
//double d_gap_density(double r, double d, double M, double alpha, double q, double rho_0){

	double qNL	= 1.04*pow(M, -3.);
	double qw	= 34.*qNL*pow(alpha*M, 0.5);
	double D	= 7.*pow(M, 1.5)/pow(alpha, 0.25);
	double A, dA, denom, q_denom, del_q, ddel_q, qr, dqr, drE;
	
	q_denom	= 1 + D*D*D*pow( pow( r/d, 1./6.) - 1., 6.);
	qr	= q*pow(q_denom, -1./3.);
	dqr	= (-q/3)*pow(q_denom, -4./3.)*D*D*D*pow( pow(r/a, 1./6.) - 1. , 5.)*pow(r/d, -5./6.)/d;
	
	if (qr < qNL){
		del_q = 1.;
		ddel_q = 0.;
	} else {
		del_q	= pow(qr/qNL, -0.5) + pow(qr/qw, 3.);
		ddel_q	= -0.5*pow(qr/qNL, -1.5)*(dqr/qNL) + 3.*pow(qr/qw, 2.)*(dqr/qw);
	}
	
	A	= (0.45/(3.*M_PI))*qr*qr*pow(M, 5.)/alpha;
	dA	= (0.45/(3.*M_PI))*2.*qr*dqr*pow(M, 5.)/alpha;
	
	denom	= 1 + A*del_q;
	drE	= -rho_0*pow(denom, -2)*(dA*del_q + A*ddel_q);
	double E = rho_0/denom;
	//return drE;
	struct Gap_Vals vals = {E, drE};
	return vals;
}

void initial( double * prim , double * x ){

   double r		= x[0];
   double phi		= x[1];
   double rx		= r*cos(phi);
   double ry		= r*sin(phi);
   double R		= sqrt( (a-rx)*(a-rx) + ry*ry );

   double mu		= q_planet/(1.+q_planet);
   /*
   double omega		= sqrt(1./(R*R*R));

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
   
   double rho_0		= 1.0;
   double nom		= f_r*K/(3.0*M_PI);
   double denom		= 1.0 + f_0*K/(3.0*M_PI);
   double rho		= rho_0 * (1.0 - (nom/denom)*sqrt(a/R));
   if(q_planet <= 0.0)
       rho = rho_0;
   //double Pp		= rho/(gam*Mach*Mach);
   //double Pp		= rho/(gam*mach_csd*mach_csd);
   */
   
   double rho_0	= 1.0;
   struct Gap_Vals vals = gap_density(R, a, mach_csd, alpha_csd, q_planet, rho_0);
   double rho = vals.rho;
   double d_rho = vals.d_r_rho;
   //double rho = gap_density(R, a, mach_csd, alpha_csd, q_planet, rho_0);
   //double d_rho = d_gap_density(R, a, mach_csd, alpha_csd, q_planet, rho_0);
   double Pp = rho * get_cs2(x, rho) / gam;

   double omega2 = (1./(R*R*R)) - (1/(R*gam*rho))*get_cs2(x, rho)*d_rho;

   double X = 0.0; 
   if( r*cos(x[1]) > 0.0 ) X = 1.0; 
   
   //double theta		= atan2(ry, a-rx);
   //double xi		= M_PI - theta - phi;
   double dr_dPhi	= -a*sin(phi);
   double dphi_dPhi	= 1.0 - (a/r)*cos(phi);

   double vr_global	= dr_dPhi*(sqrt(omega2) - sqrt(1./(a*a*a)));
   double vp_global	= dphi_dPhi*(sqrt(omega2) - sqrt(1./(a*a*a))) + rot_om;
   double vp_local	= sqrt(mu/(r*r*r));
   
   //Fixing angular velocity profile 
   double r_in		= 0.005;
   double r_out		= 0.05;
   double vr_f, vp_f, slope_p, slope_r, slope_rho;
   double rho_in	= 0.5*rho_0;
  
    double X2;

   if (r < r_in)
   {
	vp_f	= vp_local;
	vr_f	= 0;
	rho	= rho_in;
    X2 = 1.0;
   }
   else if (r_in < r && r < r_out)
   {
	slope_p	= (vp_global - vp_local)/(r_out - r_in);
	vp_f	= vp_local + slope_p*(r-r_in);

	slope_r	= (vr_global - 0)/(r_out - r_in);
	vr_f	= 0 + slope_r*(r-r_in);
	
	slope_rho = (rho - rho_in)/(r_out - r_in);
	rho	= rho_in + slope_rho*(r-r_in);
    X2 = (r_out - r) / (r_out - r_in);
   }else{
	vp_f	= vp_global;
	vr_f	= vr_global;
	rho	= rho;
    X2 = 0.0;
   }

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr_f; //-sin(xi) * (sqrt(1./(R)) - sqrt(1./(a)); //-1.5*nu/r;
   prim[UPP] = vp_f; //cos(xi) R * sqrt(1./(a*a*a));
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;
   if( NUM_N>1 ) prim[NUM_C+1] = X2;

}
