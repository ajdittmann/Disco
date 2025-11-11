#include "../paul.h"
#include "../Calc/calc.h"

static double M = 0.0;
static double gam = 0.0;
static int cyl = 0;

void setICparams( struct domain * theDomain ){
   gam = theDomain->theParList.Adiabatic_Index;
   M = theDomain->theParList.metricPar2;
   cyl = 0;
}

void initial( double * prim , double * x ){

    double r = x[0];
    double z = x[2];
    double R = sqrt(r*r+z*z);
    double R0 = 10.0;
    double rs = 10.0;
    double Mdot = 1.0;
    double b0 = 0.0;
    
    double us2 = M / (2*rs);
    double as2 = us2 / (1 - 3*us2);
    double us = -sqrt(us2);

    double rhos = -Mdot / (4*M_PI*rs*rs*us);
    double K = as2 / (gam * pow(rhos,gam-1) * (1-as2/(gam-1)));

    double a02 = (gam-1) * (1 - sqrt(1+3*as2) * fabs(1-as2/(gam-1)));
    double rho0 = pow(a02 / (gam * K * (1-a02/(gam-1))), 1.0/(gam-1));
    double h0 = 1.0 + gam/(gam-1) * K * pow(rho0, gam-1);

    //printf("M=%.6lg rs=%.6lg as=%.6lg us=%.6lg rhos=%.6lg K=%.6lg a02=%.6lg rho0=%.6lg h0=%.6lg\n",
    //        M, rs, sqrt(as2), us, rhos, K, a02, rho0, h0);

    double uRsc = bondi_rel_solve(Mdot, M, R, gam, h0, K, sqrt(as2));
    double rhoB = -Mdot / (4*M_PI*R*R*uRsc);
    double PB = K * pow(rhoB, gam);


    double u0sc = sqrt((1+uRsc*uRsc/(1-2*M/R))/(1-2*M/R));

    double uRks = uRsc;
    double u0ks = u0sc + uRsc/(R/(2*M)-1);
    double lRks = 2*M/R*u0ks + (1+2*M/R)*uRks;

    double chi = R>R0 ? 1.0 : 0.5*(1+tanh(tan(M_PI*(R/R0-0.5))));

    double rho = rhoB;
    double P = PB;
    double lR = lRks;

    if(R < rs)
    {
        rho = rhos;
        P = K * pow(rhos, gam);
        lR = -sqrt(2*M/R) / (1+sqrt(2*M/R));
    }

    if(cyl)
    {
        prim[RHO] = rho;
        prim[PPP] = P;
        prim[URR] = r/R * lR * chi;
        prim[UPP] = 0.0;
        prim[UZZ] = z/R * lR * chi;

        if(NUM_C > 5)
        {
            prim[BRR] = r/(R*R*R)*b0 * chi;
            prim[BPP] = 0.0;
            prim[BZZ] = z/(R*R*R)*b0 * chi;
        }
    }
    else
    {
        prim[RHO] = rho;
        prim[PPP] = P;
        prim[URR] = lR * chi;
        prim[UPP] = 0.0;
        prim[UZZ] = 0.0;

        if(NUM_C > 5)
        {
            prim[BRR] = b0/(R*R) * chi;
            prim[BPP] = 0.0;
            prim[BZZ] = 0.0;
        }
    }

    int q;
    for(q=NUM_C; q < NUM_Q; q++)
    {
        prim[q] = 0.0;
    }
}
