#include "../disco.h"

static int visc_flag = 0;
static double nu = 0.0;
static double r0 = 0.0;
static double dr = 0.0;
static double dr2 = 0.0;
static int prof = 0;
static double sig0 = 0.0;
static double mach = 0.0;
static double gam = 0.0;
static double M = 0.0;
static double a = 0.0;
static double a_off = 0.0;
static double q = 0.0;
static double rg = 0.0;
static double eps = 0.0;

void setICparams( struct domain * theDomain ){
    visc_flag = theDomain->theParList.visc_flag;
    nu = theDomain->theParList.viscosity;
    prof = theDomain->theParList.initPar0;
    r0 = theDomain->theParList.initPar1;
    dr = theDomain->theParList.initPar2;
    dr2 = theDomain->theParList.initPar3;
    sig0 = theDomain->theParList.initPar4;
    a_off = theDomain->theParList.initPar5;
    mach = theDomain->theParList.Disk_Mach;
    gam = theDomain->theParList.Adiabatic_Index;
    eps = theDomain->theParList.grav_eps;
    
    M = 1.0;
    q = theDomain->theParList.Mass_Ratio;
    if(strcmp(PLANETS, "cm") == 0)
        a = 1.0;
    else
        a = 0.0;
    if(strcmp(PLANETS, "bh") == 0)
        rg = M;
    else
        rg = 0.0;
}

double poly_step(double x, double a, double b, double ya, double yb)
{
    double y;
    if(x <= a)
        y = ya;
    else if(x >= b)
        y = yb;
    else
    {
        double X = (2*x - (a+b)) / (b-a);
        double Y = 15.0*X/8.0 * (1.0 - 2.0*X*X/3.0 + X*X*X*X/5.0);
        y = 0.5*(yb-ya) * Y + 0.5*(ya+yb);
    }

    return y;
}

double dpoly_step(double x, double a, double b, double ya, double yb)
{
    double dy;
    if(x <= a)
        dy = 0.0;
    else if(x >= b)
        dy = 0.0;
    else
    {
        double X = (2*x - (a+b)) / (b-a);
        double dX = 2.0/(b-a);
        //double Y = 15.0*X/8.0 * (1.0 - 2.0*X*X/3.0 + X*X*X*X/5.0);
        double dY = 15.0/8.0 * (1.0 - 2.0*X*X + X*X*X*X) * dX;
        dy = 0.5*(yb-ya) * dY;
    }

    return dy;
}

void initial( double * prim , double * x ){

   double r = x[0];
   double phi = x[1];

   double rho, drhodr;
   double rm = r0-dr;
   double rp = r0+dr;

   double sig1 = 1.0;

   double xc = r * cos(phi);
   double yc = r * sin(phi);
   double X = xc - a_off;
   double Y = yc;
   double R = sqrt(X*X + Y*Y);
   double PHI = atan2(Y, X);

   
   if(prof == 1 || prof == -1) //Gaussian
   {
       rho = sig0 + (sig1-sig0)*exp(-0.5*(R-r0)*(R-r0)/(dr*dr));
       drhodr = -(sig1-sig0)*(R-r0)/(dr*dr);
   }
   else if (prof == 2 || prof == -2)  //Smoothed top hat ~ tanh(1-r^2)
   {
       rho = sig0 + (sig1-sig0) * 0.5*(1.0+tanh(-(R-rm)*(R-rp)/ (dr2*dr2)));
       double coshr = cosh((R-rm)*(R-rp)/(dr2*dr2));
       drhodr = -(sig1-sig0) * (R-0.5*(rm+rp)) / (coshr*coshr*dr2*dr2);
   }
   else if (prof == 3 || prof == -3)  //Top Hat with cosine boundaries
   {
       double rm2 = rm - dr2;
       double rp2 = rp + dr2;
       if(R > rm && R < rp)
       {
           rho = sig1;
           drhodr = 0.0;
       }
       else if(R>rm2 && R < rm)
       {
           rho = sig0 + (sig1-sig0) * 0.5*(1+cos(M_PI*(R-rm)/dr2));
           drhodr = -0.5*M_PI*(sig1-sig0)/dr2 * sin(M_PI*(R-rm)/dr2);
       }
       else if(R>rp && R < rp2)
       {
           rho = sig0 + (sig1-sig0) * 0.5*(1+cos(M_PI*(R-rp)/dr2));
           drhodr = -0.5*M_PI*(sig1-sig0)/dr2 * sin(M_PI*(R-rp)/dr2);
       }
       else
       {
           rho = sig0;
           drhodr = 0.0;
       }
   }
   else if (prof == 4 || prof == -4)  //Top Hat with polynomial boundaries
   {
       double rm2 = rm - dr2;
       double rp2 = rp + dr2;
       if(R > rm && R < rp)
       {
           rho = sig1;
           drhodr = 0.0;
       }
       else if(R>rm2 && R < rm)
       {
           rho = poly_step(R, rm2, rm, sig0, sig1);
           drhodr = dpoly_step(R, rm2, rm, sig0, sig1);
       }
       else if(R>rp && R < rp2)
       {
           rho = poly_step(R, rp, rp2, sig1, sig0);
           drhodr = dpoly_step(R, rp, rp2, sig1, sig0);
       }
       else
       {
           rho = sig0;
           drhodr = 0.0;
       }
   }
   else if (prof == 5 || prof == -5)  //Polynomial step up
   {
       double rm2 = rm - dr2;
       //double rp2 = rp + dr2;
       if(R > rm)
       {
           rho = sig1;
           drhodr = 0.0;
       }
       else if(R>rm2 && R < rm)
       {
           rho = poly_step(R, rm2, rm, sig0, sig1);
           drhodr = dpoly_step(R, rm2, rm, sig0, sig1);
       }
       else
       {
           rho = sig0;
           drhodr = 0.0;
       }
   }
   else
   {
       rho = sig1;
       drhodr = 0.0;
   }

   double cs20 = 1.0/r0 / (mach*mach);
   double Pp0 = sig1 * cs20 / gam;

   double Pp, dPdr;
   if(prof > 0)
   {
      Pp = Pp0 * pow(rho / sig1, gam);
      dPdr = (gam*Pp/rho) * drhodr;
   }
   else
   {
       Pp = Pp0;
       dPdr = 0.0;
   }

   //double qp1 = 1+q;
   //double y = a/((1+q)*R);
   //double gr = -M/((R-2*rg)*(R-2*rg)) * (1 + 3./4. * q * y*y
   //                     + 45./64. * q*(q*q-q+1) * y*y*y*y);

   double gr = -M * R / pow(R*R + eps*eps, 1.5);
   double OM = sqrt((dPdr/rho - gr) / R);
   double VR = visc_flag ? -1.5*nu/rho/R : 0.0;
   if(R < 6*rg)
   {
       OM = 1.5*sqrt(M*6*rg)/(R*R);
   }

   double vx = VR * cos(PHI) - R*OM * sin(PHI);
   double vy = VR * sin(PHI) + R*OM * cos(PHI);

   double vr = vx * cos(phi) + vy * sin(phi);
   double om = (-vx * sin(phi) + vy * cos(phi)) / r;

   double Q = 0.0; 
   if( R>rm && R<rp ) Q = 1.0; 
   double Q2 = 0.0; 
   if( cos(PHI) > 0 && Q > 0) Q2 = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr;
   prim[UPP] = om;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = Q;
   if( NUM_N>1 ) prim[NUM_C+1] = Q2;
}
