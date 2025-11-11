#include "../paul.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;
static double M = 0.0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
   M = theDomain->theParList.metricPar2;
}

void initial(double *prim, double *x)
{
    double r = x[0];
    double phi = x[1];
    double z = x[2];
    double R = sqrt(r*r+z*z);
    double st = r/R;
    double ct = z/R;

    double g00 = -1 + 2*M/R;
    double g0R = 2*M/R;
    double gRR = 1 + 2*M/R;
    double gtt = R*R;
    double gpp = r*r;

    double u0eq, uReq, upeq, u0pol, uRpol, uppol;

    if(R > 6*M)
    {
        u0eq = sqrt(R/(R-3*M));
        uReq = 0.0;
        upeq = sqrt(M/(R*R*R - 3*M*R*R));
    }
    else
    {
        double x = 6*M/R - 1.0;
        u0eq = 2.0 * (sqrt(2.0)*R - M*sqrt(x*x*x)) / (3.0*(R-2.0*M));
        uReq = -sqrt(x*x*x) / 3.0;
        upeq = 2.0*sqrt(3.0)*M/(R*R);
    }
    u0pol = (R*R - 2*M*sqrt(2.0*M*R)) / (R*R - 2*M*R);
    uRpol = -sqrt(2*M/R);
    uppol = 0.0;

    double lReq, lpeq, lRpol, lppol;

    lReq = g0R*u0eq + gRR*uReq;
    lpeq = gpp*upeq;
    lRpol = g0R*u0pol + gRR*uRpol;
    lppol = gpp*uppol;


    double lR = ct*ct*lRpol + st*st*lReq;
    double lp = ct*ct*lppol + st*st*lpeq;
    double lt = 0.0;

    double lr = st * lR - ct/R * lt;
    double lz = ct * lR + st/R * lt;

    //lr = 0.0;

    double uIsco2 = 1.0/3.0;
    double u_cs2 = uIsco2 / (Mach*Mach);
    double cs2 = u_cs2 / (1.0 + u_cs2);

    double rho = 1.0;
    double Pp = rho / (1.0/cs2 - gam/(gam-1));

    prim[RHO] = rho;
    prim[PPP] = Pp;
    prim[URR] = lr;
    prim[UPP] = lp;
    prim[UZZ] = lz;

    if(NUM_C > 5)
    {
        int q;
        for(q=5; q<NUM_C; q++)
            prim[q] = 0.0;
    }
    
    if( NUM_N>0 ) 
        prim[NUM_C] = r>10.0 ? 1.0 : 0.0;
    if(NUM_N > 1)
        prim[NUM_C+1] = r*cos(phi)>0 ? 1.0 : 0.0;
    if(NUM_N > 2)
    {
        double x0 = 10.0;
        double y0 = 0.0;
        double R0 = 1.0;
        double dx = r*cos(phi) - x0;
        double dy = r*sin(phi) - y0;
        prim[NUM_C+2] = sqrt(dx*dx+dy*dy)<R0 ? 1.0 : 0.0;
    }
}
