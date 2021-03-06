#include "../paul.h"

static double gam = 0.0;
static double visc = 0.0;
static int alpha_flag = 0;
static struct planet *thePlanets = NULL;
static double Mach = 0.0;
static int Npl = 0;
static double massq = 0.0;
static double xi = 20.0;
static double rin = 0.0;
static double redge = 0.0;
static double rswitch = 0.0;
static double epsfl = 0.0;
static int cs2choice = 0;

static double a = 1.0;
static double om_bary = 0.0;

double get_cs2(double *);
double phigrav( double , double , double , int); //int here is type
double fgrav( double , double , double , int);

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    visc = theDomain->theParList.viscosity;
    alpha_flag = theDomain->theParList.alpha_flag;
    Mach = theDomain->theParList.Disk_Mach;
    massq = theDomain->theParList.Mass_Ratio;
    thePlanets = theDomain->thePlanets;
    Npl = theDomain->Npl;
    epsfl = theDomain->theParList.Density_Floor;

    //a = theDomain->theParList.RotD; //only true for  q->0
    om_bary = theDomain->theParList.RotOmega;


    //xi = theDomain->theParList.initPar1;
    //rin = theDomain->theParList.initPar2;
    rin = 2.0*theDomain->theParList.initPar2;
    redge = theDomain->theParList.initPar3;

    cs2choice = theDomain->theParList.Cs2_Profile;


}

void initial(double *prim, double *x)
{
   double r     = x[0];
   double phi       = x[1];
   double rx        = r*cos(phi);
   double ry        = r*sin(phi);
   double R     = sqrt( (a-rx)*(a-rx) + ry*ry );

   double p1r = thePlanets[1].r;
   double p1phi = thePlanets[1].phi;
   double p1rx = p1r*cos(p1phi);
   double p1ry = p1r*sin(p1phi);
   double d1p = sqrt( (rx-p1rx)*(rx-p1rx) + (ry-p1ry)*(ry-p1ry));

   double p2r = thePlanets[2].r;
   double p2phi = thePlanets[2].phi;
   double p2rx = p2r*cos(p2phi);
   double p2ry = p2r*sin(p2phi);
   double d2p = sqrt( (rx-p2rx)*(rx-p2rx) + (ry-p2ry)*(ry-p2ry));

   double homtot = 0.0;
   homtot +=  fgrav( thePlanets[0].M, R , thePlanets[0].eps, thePlanets[0].type)/R;
   double cs2 = get_cs2(x);

   double alpha = visc;
   if (alpha_flag == 0) alpha = alpha*sqrt(homtot)/cs2;

   double dphi = 0.0;
   if (cs2choice == 5){
     //should probably double check this before using...
     dphi -= phigrav( thePlanets[0].M , R , thePlanets[0].eps , thePlanets[0].type )/R;
   }
   double nu = visc;
   if (alpha_flag == 0) nu = alpha*cs2/homtot;
   double epsfact = exp(-pow(r/rin, xi));

   double Sigma = epsfact+epsfl;
   double dSigma = -xi*exp(-pow(r/rin,xi))*pow(r/rin, xi-1.0)/rin;

   //double omega_glob = sqrt(fabs(1.0/(R*R*R) + (cs2*dSigma/gam + Sigma*dphi/gam)/(R*Sigma)));

   //double vr_glob = -1.5*nu/R;
   //double vphi_glob = R*omega_glob;

   //Sigma = 1.0; //testing
   //vphi_glob = 0.0; //testing
   //vr_glob = 0.0; //testing
   //double efact = exp(-pow((r/redge),-xi));
   double rho = Sigma; //*efact; // + epsfl;
   //if (r < rin) rho*=0.01;
   double P = rho*cs2/gam;

   //double dr_dPhi	= -a*sin(phi);
   //double dphi_dPhi	= 1.0 - (a/r)*cos(phi);
   //double dr_dR = R/r - (R*R+a*a-r*r)/(2*R*r);
   //double dphi_dR = (r*sin(phi)/R - dr_dR*sin(phi))/(r*cos(phi));

   //double vr	= dr_dPhi*(omega_glob - sqrt(1./(a*a*a))) + dr_dR*vr_glob;
   //double vp	= dphi_dPhi*(omega_glob - sqrt(1./(a*a*a))) + om_bary + dphi_dR*vr_glob;

   double vr	= -1.5*nu/r;
   double vp = sqrt(fabs((massq/(1.0+massq))/(r*r*r) + (cs2*dSigma/gam + Sigma*dphi/gam)/(r*Sigma)))+om_bary;

   prim[RHO] = rho;
   prim[PPP] = P;
   prim[URR] = 0.0;
   prim[UPP] = vp;
   prim[UZZ] = 0.0;

   int q2;
   for(q2 = 5; q2 < NUM_Q; q2++)

       prim[q2] = 0.0;
}
