#include "../disco.h"

static double gam = 0.0;
static double g = 0.0;
static double H0 = 0.0;
static double a = 0.0;
static double L = 0.0;
static double x0 = 0.0;
static double gamInt = 0.0;
static double rhoatm = 0.0;
static double Patm = 0.0;
static double zbot = 0.0;
static int integrateVert = 0;

void get_vec_from_xyz(double *x, double *vxyz, double *v);
void get_vec_contravariant(double *x, double *v, double *vc);
    
void setICparams( struct domain * theDomain )
{
    g = 1.0;
    gam = theDomain->theParList.Adiabatic_Index;
    zbot = theDomain->theParList.zmin;
    if(theDomain->Nz == 1)
        integrateVert = 1;
    H0 = theDomain->theParList.initPar1;
    a = theDomain->theParList.initPar2;
    L = theDomain->theParList.initPar3;
    x0 = theDomain->theParList.initPar4;
    gamInt = theDomain->theParList.initPar5;
    rhoatm = theDomain->theParList.initPar6;
    Patm = theDomain->theParList.initPar7;
}

void initial( double * prim , double * xyz )
{
    // This is a bit of a weird setup. The idea is to model a box of gas
    // with a hard bottom in a uniform vertical gravitational field. BUT
    // we're using Disco's built in 2D gravity, which produces fields in the
    // cylindrical r direction.  So we're swapping x and z.  Goes without
    // saying this only makes sense in the cartesian setup.

    double x = xyz[0];
    double z = xyz[2];

    //Assuming isentropic gas, so P = k*rho^gam
    double k = 1.0;

    double sigref = pow((gamInt-1)/gamInt * pow(g/k, 1./gamInt) * H0, 
                        gamInt/(gamInt-1));

    double X = (x-x0)/L;

    double sig;
    if(fabs(X) < 1.0)
    {
        double f = 1.0-X*X;
        sig = sigref*(1 + a*f*f*f*f);
    }
    else
        sig = sigref;

    double H = gamInt/(gamInt-1) * pow(k/g, 1./gamInt)
                * pow(sig, (gamInt-1)/gamInt);

    double P = pow((gamInt-1)/gamInt * g * pow(k,-1./gamInt) * (H-(z-zbot)), 
                    gamInt/(gamInt-1));
    double rho = pow(P/k, 1./gamInt);

    if(P <= Patm || rho <= rhoatm || z-zbot >= H) 
    {
        P = Patm;
        rho = rhoatm;
    }

    if(integrateVert)
    {
        rho = sig;
        P = sig * (gamInt-1)/(2*gamInt-1) * g*H;
    }
    
    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = 0.0;
    prim[UPP] = 0.0;
    prim[UZZ] = 0.0;

    int q;
    for(q = NUM_C; q < NUM_Q; q++)
        prim[q] = 0.0;
}
