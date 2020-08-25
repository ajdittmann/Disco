#include "../paul.h"

static double gam = 0.0;
static double visc = 0.0;
static int alpha_flag = 0;
static struct planet *thePlanets = NULL;
static double Mach = 0.0;
static int Npl = 0;
static double massq = 0.0;
static double xi = 0.0;
static double rin = 0.0;
static double redge = 0.0;
static double rswitch = 0.0;

double get_cs2(double *);

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    visc = theDomain->theParList.viscosity;
    alpha_flag = theDomain->theParList.alpha_flag;
    Mach = theDomain->theParList.Disk_Mach;
    massq = theDomain->theParList.Mass_Ratio;
    thePlanets = theDomain->thePlanets;
    Npl = theDomain->Npl;
    xi = theDomain->theParList.initPar1;
    rin = theDomain->theParList.initPar2;
    redge = theDomain->theParList.initPar3;
    rswitch = theDomain->theParList.rmax;
    rswitch = rswitch*0.5;
}

void initial(double *prim, double *x)
{
    double epsfl = 0.0001;
    double r = x[0];
    double R = r + epsfl;
    double phi = x[1];

    double cs2 = get_cs2(x);

    double rho, efact, ztfact1, ztfact, fth, dfth;

    //double om = 1.0;
    double om = pow(R,-1.5);
    if (R<0.1) om = pow(0.1, -1.5);
    int np;
    double alpha = visc;
    double nu = visc;
    if (alpha_flag == 1){
      if (Npl < 2){
          nu = alpha*cs2/sqrt(om);
      }
      else{
        double omtot = 0;
        double cosp, sinp, px, py, dx, dy, gx, gy, mag;
        gx = r*cos(phi);
        gy = r*sin(phi);
        for(np = 0; np<Npl; np++){
          cosp = cos(thePlanets[np].phi);
          sinp = sin(thePlanets[np].phi);
          px = thePlanets[np].r*cosp;
          py = thePlanets[np].r*sinp;
          dx = gx-px;
          dy = gy-py;
          mag = dx*dx + dy*dy + thePlanets[np].eps*thePlanets[np].eps;
          omtot +=	thePlanets[np].M*pow(mag, -1.5);
        }  	
        nu = alpha*cs2/sqrt(omtot);
        om = sqrt(omtot);
      }
    }

    double phitot = 0.0;
    double dphitot = 0.0;
    for (np = 0; np<Npl; np++){
      double denom = thePlanets[np].eps*thePlanets[np].eps + thePlanets[np].r*thePlanets[np].r;
      denom += r*r - r*thePlanets[np].r*cos(phi - thePlanets[np].phi);
      double sqdenom = sqrt(denom);
      phitot -= thePlanets[np].M/sqdenom;
      dphitot += thePlanets[np].M*(r - thePlanets[np].r*cos(phi - thePlanets[np].phi))/(denom*sqdenom);
    }

    double sig0 = 1.0;
    efact = exp(-pow((R/redge),-xi));

    //double rho, efact, ztfact1, ztfact, fth, dfth;

    ztfact1 = 1.0 - sqrt(rin/(R));
    fth = 0.5*(1.0 + tanh(r - rswitch));
    dfth = cosh(r - rswitch);
    dfth = 0.5/(dfth*dfth);
    ztfact = ztfact1*(1.0 - fth) + fth;
    ztfact = fmax(ztfact, epsfl);

    rho = sig0*ztfact*efact + epsfl;
    double drho = efact*(0.5*sqrt(rin)/(R*sqrt(R)))*(1.0 - fth);
    drho = drho + dfth + dfth*ztfact1*(fth-1.0) + ztfact*(xi*redge*redge/(R*R*R))*efact;
 
    double v = -1.5*nu/(R*ztfact);
    double P = -rho*phitot/(Mach*Mach);

    double multom = 1.0 + 0.75*massq/(R*R*(1.0 + massq)*(1.0 + massq));
    double addom = rho*dphitot + phitot*drho;
    addom *= -1.0/(Mach*Mach*R*rho);
    om = sqrt(fabs(om*om*multom + addom));

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = v;
    prim[UPP] = om;
    prim[UZZ] = 0.0;

    int q2;
    for(q2 = 5; q2 < NUM_Q; q2++)

        prim[q2] = 0.0;
}