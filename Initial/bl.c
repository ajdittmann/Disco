#include "../disco.h"
#include "../omega.h"

static double M = 0.0;
static double as = 0.0;
static double gam = 0.0;
static double r0 = 0.0;
static double dr = 0.0;
static double sig0 = 0.0;
static double sig_atm = 0.0;
static double mach = 0.0;
static int isothermal = 0;

double x_func_sin2(double x, void *args);
double r_func_schw_sc(double r, void *args);
double int_sim_ad(double a, double b, double tol, 
                    double (*f)(double, void *), void *args);

void setICparams( struct domain * theDomain ){
   gam = theDomain->theParList.Adiabatic_Index;
   M = theDomain->theParList.metricPar2;
   as = theDomain->theParList.metricPar3;
   r0 = theDomain->theParList.initPar1;
   dr = theDomain->theParList.initPar2;
   sig0 = theDomain->theParList.initPar3;
   sig_atm = theDomain->theParList.initPar4;
   mach = theDomain->theParList.Disk_Mach;
   isothermal = theDomain->theParList.isothermal_flag;
}

void initial( double * prim , double * x ){

   double rnd1 = (double)rand()/(double)(RAND_MAX);
   double rnd2 = (double)rand()/(double)(RAND_MAX);

   double r   = x[0];
   double phi = x[1];
   double z   = 0.0; //x[2];

   double rp = r0 + 0.5*dr;
   double rm = r0 - 0.5*dr;

   double vp, u0, sig, pi;

   double cs20;
   if(isothermal)
      cs20 = get_cs2(x);
   else
      cs20 = M/(r0 * mach*mach);

   double pi0 = sig0 * cs20/gam;

   if(r > rp)
   {
       vp = sqrt(M/(r*r*r));
       u0 = 1.0/sqrt(1-3*M/r);
       sig = sig0;
       pi = pi0;
   }
   else if(r > rm)
   {
       vp = (r-rm)/dr * sqrt(M/(rp*rp*rp));
       u0 = 1.0/sqrt(1-2*M/r-r*r*vp*vp);
       double args[4] = {M, rm, rp, 0.0};
       double Ir = int_sim_ad(r, rp, 1.0e-8, &r_func_schw_sc, args);
       if(isothermal)
       {
           sig = sig0 * exp(-Ir/cs20);
           pi = sig*cs20/gam;
       }
       else
       {
           sig = sig0 * pow(1.0 - (gam-1)*Ir/cs20, 1.0/(gam-1));
           pi = pi0 * pow(sig/sig0, gam);
       }
   }
   else
   {
       vp = 0.0;
       u0 = 1.0/sqrt(1-2*M/r);
       double args[4] = {M, rm, rp, 0.0};
       double Ir = int_sim_ad(r, rp, 1.0e-8, &r_func_schw_sc, args);
       if(isothermal)
       {
           sig = sig0 * exp(-Ir/cs20);
           pi = sig*cs20/gam;
       }
       else
       {
           sig = sig0 * pow(1.0 - (gam-1)*Ir/cs20, 1.0/(gam-1));
           pi = pi0 * pow(sig/sig0, gam);
       }
   }

   if (sig < sig_atm)
       sig = sig_atm;

   //vr = (1e-3*rnd1 - 5e-4)

   prim[RHO] = sig;
   prim[PPP] = pi;
   prim[URR] = 0.0;
   prim[UPP] = r*r*u0*vp;
   prim[UZZ] = 0.0;

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}

double r_func_schw_sc(double r, void *args)
{
    double M = ((double *)args)[0];
    double rm = ((double *)args)[1];
    double rp = ((double *)args)[2];
    double om = ((double *)args)[3];

    double vp = (r-rm)/(rp-rm) * (sqrt(M/(rp*rp*rp)) - om) + om;
    if(r < rm)
        vp = om;

    return (M - r*r*r*vp*vp) / (2*M*r - r*r + r*r*r*r*vp*vp);
}

double x_func_sin2(double x, void *args)
{
    return sin(x)*sin(x);
}

double int_sim_ad(double a, double b, double tol, double (*f)(double, void *), void *args)
{
    int N1 = 4;

    double S = 0.0;
    double T = 0.0;
    double I = 0.0;
    double I1 = 0.0;

    int i, n;

    int N = N1;
    double h = (b-a) / N;

    S = f(a, args) + f(b, args);
    for(i=2; i<N; i+=2)
        S += 2*f(a+i*h, args);
    for(i=1; i<N; i+=2)
        T += f(a+i*h, args);
    S /= 3.0;
    T *= 2.0/3.0;
    I1 = h*(S + 2*T);

    double err = 1.0e10 * tol;

    while(fabs(err) > tol)
    {
        N *= 2;
        h *= 0.5;
        I = I1;

        S += T;
        T = 0.0;
        for(i=1; i<N; i+=2)
            T += f(a+i*h, args);
        T *= 2.0/3.0;
        I1 = h*(S + 2*T);

        err = (I1-I)/15.0;
    }

    return I1;
}
