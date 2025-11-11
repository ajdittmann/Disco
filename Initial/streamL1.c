#include "../paul.h"

// Taken directly from Lubow & Shu 1975

static double gam = 0.0;
static double eps = 0.0;
static double Mdot = 0.0;
static double R0 = 0.0;
static double fac_atm = 0.0;
static double cs_atm = 0.0;
static double q = 0.0;
static double Omega0 = 0.0;
static int twoD = 0;
static double xL1 = 0.0;

void get_xyz(double *x, double *xyz);
void get_rpz(double *x, double *xyz);
void get_vec_from_xyz(double *x, double *vxyz, double *v);
void get_vec_rpz(double *x, double *v, double *vrpz);
void get_vec_from_rpz(double *x, double *vrpz, double *v);
void get_vec_contravariant(double *x, double *v, double *vc);

double calc_L1(double Q)
{
    //Calculate the distance from the center of mass to L1, where 
    //positive denotes towards the primary (M1) and negative denotes
    //towards the secondary (M2). In DISCO M1 typically begins on the
    //+'ve x-axis (phi=0) and M2 on the -'ve x-axis (phi=pi).
    //q = M2/M1
    //
    //Assume the total mass G*M and the binary separatation a are 1.

    double x = 0.0;
    double M1 = 1.0/(1+Q);
    double M2 = Q/(1+Q);
    double a1 = Q/(1+Q);
    double a2 = 1.0/(1+Q);
    double xa = -a2;
    double xb = a1;

    double dx = 10.0;
    double f, df;

    int i = 0;
    while(fabs(dx) > 1.0e-12 && i<100)
    {
        double d1 = x-a1;
        double d2 = x+a2;
        f = M1/(d1*d1) - M2/(d2*d2) + x;
        df = -2*M1/(d1*d1*d1) + 2*M2/(d2*d2*d2) + 1;

        if(f > 0.0)
            xb = x;
        if(f < 0.0)
            xa = x;

        dx = -f/df;

        //printf("%02d: x=%.6lg f=%.6lg df=%.6lg\n", i, x, f, df);
        //printf("      dx=%.6lg xa=%.6lg xb=%.6lg (xb-xa)=%.6lg\n", 
        //            dx, xa, xb, xb-xa);

        if(x+dx <= xa || x+dx >= xb)
        {
            dx = 0.5*(xa+xb)-x;
            x = 0.5*(xa+xb);
        }
        else
            x += dx;
        i += 1;
    }

    printf("From M2: CM=%.6lf L1=%.6lf\n",a2, a2+x); 

    return x;
}
    
void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    eps = theDomain->theParList.initPar1;
    Mdot = theDomain->theParList.initPar2;
    R0 = theDomain->theParList.initPar3;
    fac_atm = theDomain->theParList.initPar4;
    cs_atm = theDomain->theParList.initPar5;
    q = theDomain->theParList.Mass_Ratio;
    Omega0 = theDomain->theParList.RotOmega;
    if(theDomain->Nz == 1)
        twoD = 1;
    else
        twoD = 0;
    xL1 = calc_L1(q);
}

void initial( double * prim , double * X )
{

    //Assume origin is on M2 and M1 is located on +'ve x axis.
    //
    // Our M2 is L&S's star D (detached - has disk)
    // Our M1 is L&S's star C (contact - donor)
    
    double xyz[3];
    get_xyz(X, xyz);

    //x,y,z are equivalent to Lubow & Shu: origin at L1, +'ve x-axis
    //towards star D.
    double x = -(xyz[0] - 1./(1+q) - xL1);
    double y = -xyz[1];
    double z = xyz[2];

    double r1 = sqrt(x*x + y*y) / eps;
    double z1 = z/eps;
    double th = atan2(y, x);
   
    double mu = q/(1+q);  // L&S: mu = MD/M = (our) M2/M
    double XL1 = -xL1;    // L&S: +'ve XL1 is towards MD (our M2)
    double f1 = fabs(XL1-1+mu);
    double f2 = fabs(XL1+mu);
    double A = mu/(f1*f1*f1) + (1-mu)/(f2*f2*f2);  //eq 13

    double thB = 0.5*acos(-(2+A)/(3*A)); //eq 16

    double iA = 1.0/A;
    
    double sin2thS = -sqrt((8.*iA/9.)*(1-2*iA+3*sqrt(1-8*iA/9.)));//24
    double thS = 0.5*asin(sin2thS);
    double US = -0.75*A*sin2thS;  //eq 23b

    double TH = 0.5*(US*US-A+2)*(th-thS)*(th-thS);  // eq 25a
    double U = US + 2*(th-thS);  // eq 25b
    double V = -US*(th-thS);  // eq 25c
    double CS = sqrt(A*(US*US-A+2)) / (2*M_PI*US);  // eq 26

    //printf("q=%.4lf XD=%.4lf XC=%.4lf XL1=%.4lf A=%.4lf\n", q, 1-mu, -mu,
    //        XL1, A);
    //printf("   ths=%.2lf thB=%.2lf US=%.3lf CS=%.4lf\n",
    //            thS*180.0/M_PI, thB*180.0/M_PI, US, CS);

    // Atmosphere
    double rho0;
    if(twoD)
        rho0 = CS * sqrt(2*M_PI/A) / (eps*eps);
    else
        rho0 = CS / (eps*eps*eps);

    double rho, P, v[3], q1;

    if(fabs(th) < thB && r1*eps < R0)
    {
        // In the stream!

        if(twoD)
        {
            double rho2;
            rho2 = CS/r1 * sqrt(2*M_PI/A) * exp(-r1*r1*TH); //eq 15,20,21
            rho = rho2 / (eps*eps) + rho0*fac_atm;
        }
        else
        {
            double rho3;
            rho3 = CS/r1 * exp(-r1*r1*TH-0.5*A*z1*z1); //eq 15,20,21
            rho = rho3/(eps*eps*eps) + rho0*fac_atm;
        }
        P = rho*eps*eps;

        double ur = eps * r1*U; //eq 19
        double uth = eps * r1*V; //eq 19

        double uXYZ[3] = {ur*cos(th)-uth*sin(th), ur*sin(th)+uth*cos(th), 0.0};
        // Put back in Disco Orientation.
        uXYZ[0] = -uXYZ[0];
        uXYZ[1] = -uXYZ[1];
        get_vec_from_xyz(X, uXYZ, v);
        //adjust for rotating frame
        double rpz[3], vrpz[3];
        get_rpz(X, rpz);
        get_vec_rpz(X, v, vrpz);
        vrpz[1] += rpz[0] * Omega0;
        get_vec_from_rpz(X, vrpz, v);

        // v is now velocity for Disco prims
        get_vec_contravariant(X, v, v);

        q1 = 1.0;
    }
    else
    {
        rho = fac_atm * rho0;
        P = rho * cs_atm*cs_atm / gam;

        double vxyz[3] = {0.0, 0.0, 0.0};
        get_vec_from_xyz(X, vxyz, v);
        get_vec_contravariant(X, v, v);

        q1 = 0.0;
    }

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = v[0];
    prim[UPP] = v[1];
    prim[UZZ] = v[2];

    if(NUM_N > 0)
    {
       double s = log(P * pow(rho, -gam)) / (gam - 1.0);
       prim[NUM_C] = s;
        
    }

    int q;
    for(q = NUM_C+1; q < NUM_Q; q++)
    {
        if(q == NUM_C+1)
            prim[q] = q1;
        else
            prim[q] = 0.0;
    }
}
