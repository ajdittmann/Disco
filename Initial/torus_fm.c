// Fishbone & Moncrief 1976 Astrophysical Journal 207, 962-976
// http://adsabs.harvard.edu/abs/1976ApJ...207..962F

#include "../disco.h"
#include "../Hydro/frame.h"

static double M = 0.0;
static double a = 0.0;
static double gam = 0.0;
static double Rin = 0.0;
static double Rmax = 0.0;
static double rho_atm = 0.0;
static double B0 = 0.0;
static int field_choice = 0;

void get_rpz(double *x, double *rpz);
void get_vec_from_rpz(double *x, double *vrpz, double *v);
void get_vec_contravariant(double *x, double *v, double *vc);

void setICparams( struct domain * theDomain ){
    gam = theDomain->theParList.Adiabatic_Index;
    //M = theDomain->theParList.metricPar2;
    M = 1.0;
    a = theDomain->theParList.metricPar3;
    field_choice = theDomain->theParList.initPar0;
    Rin = theDomain->theParList.initPar1;
    Rmax = theDomain->theParList.initPar2;
    rho_atm = theDomain->theParList.initPar3;
    B0 = theDomain->theParList.initPar4;
}

void initial( double * prim , double * x ){

    double rpz[3];
    get_rpz(x, rpz);
    double r   = rpz[0];
    double z  = rpz[2];
    double R = sqrt(r*r + z*z);

    double k;
    if(Rmax > 2*Rin)
        Rmax = 2*Rin;
    if(Rmax < Rin)
        Rmax = Rin;
    k = Rmax/Rin;

    //double l2 = M*Rmax;
    //double h = M/R - 0.5*l2/(r*r) - 0.5*M*(2-k)/Rin;
    //double hmax = 0.5*M*(k-1)*(k-1)/Rmax;
    double rs = 2*M;
    double l2 = M*Rmax*Rmax*Rmax/((Rmax-rs)*(Rmax-rs));
    double h = M/(R-rs) - 0.5*l2/(r*r) - M/(Rin-rs) + 0.5*l2/(Rin*Rin);
    double hmax = M/(Rmax-rs) - 0.5*l2/(Rmax*Rmax)
                    - M/(Rin-rs) + 0.5*l2/(Rin*Rin);
    double rho = pow(h/hmax, 1.0/(gam-1.0));
    double l = sqrt(l2);
    double om = l/(r*r);
    double q = 1.0;

    if(h < 0.0 || rho < rho_atm || R < Rin)
    {
        //h = M/R;
        //rho = rho_atm * pow(Rin/R, 1.0/(gam-1.0));
        h = M/(R-rs);
        rho = rho_atm * pow((Rin-rs)/(R-rs), 1.0/(gam-1.0));
        q = 0.0;
        om = 0.0;
    }

    double P = (gam-1)/gam * rho*h;

    double Br, Bp, Bz;

    if(NUM_C > BZZ)
    {
        if(field_choice == 0)
        {
            Br = 0.0;
            Bp = 0.0;
            Bz = B0;
        }
        else if(field_choice == 1)
        {
            // A_phi ~ max(rho/rho_max - 0.2, 0)
            // B^r   = -d_z A_phi
            // B^z   =  d_r A_phi
            // B^phi =          0
            // smooth max() function as (0.2^alpha + rho^alpha)^(1/alpha)
            
            double alpha = 4;
            double rho_cut = 0.2;
            double dfrho = pow(1 + pow(rho_cut/rho, alpha), 1./alpha - 1);
            if (h < 0.0 || rho<rho_atm || R<Rin)
                dfrho = 0.0;
            //double dhdr = -M*r/(R*R*R) + M*Rmax/(r*r*r);
            //double dhdz = -M*z/(R*R*R);
            dfrho = rho>0.2 ? 1.0 : 0.0;
            double dhdr = -M*r/((R-rs)*(R-rs)*R) + l2/(r*r*r);
            double dhdz = -M*z/((R-rs)*(R-rs)*R);
            double drhodh = 1.0/(gam-1.0) * rho/h;
            double dAdr = B0*Rmax*Rmax * dfrho*drhodh*dhdr;
            double dAdz = B0*Rmax*Rmax * dfrho*drhodh*dhdz;
           
            Br = -dAdz / r;
            Bp =  0.0;
            Bz =  dAdr / r;
        }
        else
        {
            Br = 0.0;
            Bp = 0.0;
            Bz = 0.0;
        }
    }

    double v[3], u[3], B[3];
    double vrpz[3] = {0.0, r*om, 0.0};
    double Brpz[3] = {Br, 0.0, Bz};
    get_vec_from_rpz(x, vrpz, v);
    get_vec_contravariant(x, v, u);
    get_vec_from_rpz(x, Brpz, B);

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = u[0];
    prim[UPP] = u[1];
    prim[UZZ] = u[2];

    if(NUM_C > BZZ)
    {
        prim[BRR] = B[0];
        prim[BPP] = B[1];
        prim[BZZ] = B[2];
    }

    if(NUM_N > 0)
        prim[NUM_C] = q;
}

