#include "../disco.h"

static int rel = 0;
static double gam = 0.0;
static double rho_ref = 0.0;
static double P_ref = 0.0;
static double v_ref = 0.0;
static double L = 0.0;
static double a = 0.0;
static double x0 = 0.0;
static double theta = 0.0;

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    P_ref = theDomain->theParList.initPar1; //100 in RAM
    a = theDomain->theParList.initPar2;     //1.0 in RAM
    x0 = theDomain->theParList.initPar3;  
    L = theDomain->theParList.initPar4;     //0.3 in RAM
    rho_ref = 1.0;
    v_ref = 0.0;
    if(strcmp(HYDRO, "greuler") == 0 || strcmp(HYDRO, "grmhd") == 0
            || strcmp(HYDRO, "grmhd2") == 0)
        rel = 1;
    else
        rel = 0;
}

void initial(double *prim, double *x)
{
    double r = x[0];
    double phi = x[1];

    double X = (r*cos(phi) - x0) / L;

    double rho, v, ur, up, P;

    double f = X*X-1.0;
    if(fabs(X) < 1.0)
    {
        rho = rho_ref*(1.0 + a*f*f*f*f);
        P = P_ref*pow(rho/rho_ref, gam);
    }
    else
    {
        rho = rho_ref;
        P = P_ref;
    }


    if(rel)
    {
        double v, w;

        if(fabs(X) < 1.0)
        {
            double cs, cs_ref, Jm, js, sgmo;
            cs = sqrt(gam*P/(rho + gam/(gam-1)*P));
            cs_ref = sqrt(gam*P_ref/(rho_ref+gam/(gam-1)*P));

            sgmo = sqrt(gam-1.0);

            Jm = 0.5 * log((1.0+v_ref)/(1.0-v_ref))
                - log((sgmo+cs_ref)/(sgmo-cs_ref)) / sgmo;
            js = exp(2*Jm + 2*log((sgmo+cs)/(sgmo-cs))/sgmo);

            v = (js - 1.0) / (js + 1.0);
        }
        else
            v = v_ref;

        w = 1.0 / sqrt(1-v*v);

        ur = w*v * cos(phi);
        up = -w*v * sin(phi) * r;
    }
    else
    {
        double v;
        if(fabs(X) < 1.0)
        {
            double cs, cs_ref;
            cs = sqrt(gam*P/rho);
            cs_ref = sqrt(gam*P_ref/rho_ref);

            v = v_ref - 2*(cs_ref - cs) / (gam - 1);
        }
        else
            v = v_ref;
        ur = v * cos(phi);
        up = -v * sin(phi) / r;
    }

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = ur;
    prim[UPP] = up;
    prim[UZZ] = 0.0;

    int q;
    for(q = 5; q < NUM_Q; q++)
        prim[q] = 0.0;
}
