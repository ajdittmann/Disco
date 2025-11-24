#include "../disco.h"
#include "../geometry.h"
#include "../planet.h"

#define N_AUX_PER_PLANET 13

static double gamma_law = 0.0;
static int Npl = 0;

/*
 * 
 */

void setReportParams(struct domain *theDomain)
{
    gamma_law = theDomain->theParList.Adiabatic_Index;
    Npl = theDomain->Npl;

    if(Npl < 2)
    {
        printf("This report needs AT LEAST 2 planets to work!\n");
        exit(1);
    }
}

int num_shared_reports()
{
    return Npl*5;
}

int num_distributed_aux_reports()
{
    return Npl * N_AUX_PER_PLANET;
}

int num_distributed_integral_reports()
{
    return 9;
}

void get_shared_reports(double *Q, struct domain *theDomain)
{
    for (int k=0; k<Npl; k++){
        Q[k*5] = theDomain->thePlanets[k].M;
        Q[k*5+1] = theDomain->thePlanets[k].r;
        Q[k*5+2] = theDomain->thePlanets[k].phi;
        Q[k*5+3] = theDomain->thePlanets[k].vr;
        Q[k*5+4] = theDomain->thePlanets[k].omega;
    }
}

void get_distributed_aux_reports(double *Q, struct domain *theDomain)
{
    double *pl_aux1 = theDomain->pl_aux;
    double *pl_aux2 = theDomain->pl_aux + NUM_PL_AUX;
    double *pl_aux3 = theDomain->pl_aux + NUM_PL_AUX*2;

    int idx_aux[N_AUX_PER_PLANET] = {
                PL_SNK_M, PL_GRV_JZ, PL_SNK_JZ, PL_GRV_PX, PL_GRV_PY,
                PL_SNK_PX, PL_SNK_PY, PL_GRV_K, PL_SNK_K,
                PL_SNK_MX, PL_SNK_MY, PL_SNK_SZ, PL_SNK_U};

    int q;
    for(q=0; q<N_AUX_PER_PLANET; q++)
    {
        Q[Npl*q + 0] = pl_aux1[idx_aux[q]];
        Q[Npl*q + 1] = pl_aux2[idx_aux[q]];
        Q[Npl*q + 2] = pl_aux3[idx_aux[q]];
    }
}

void get_distributed_integral_reports(const double *x, const double *prim,
                                      double *Q, struct domain *theDomain)
{
    double rpz[3];
    get_rpz(x, rpz);

    double r = rpz[0];
    double phi = rpz[1];
    double z = rpz[2];

    double rho = prim[RHO];

    double V[3] = {prim[URR], prim[UPP], prim[UZZ]};
    get_vec_covariant(x, V, V);
    double Vrpz[3];
    get_vec_rpz(x, V, Vrpz);

    double cosp = cos(phi);
    double sinp = sin(phi);
    double cos2p = (cosp-sinp)*(cosp+sinp);
    double sin2p = 2*cosp*sinp;

    double xyz[3] = {r*cosp, r*sinp, z};

    double Fxyz[3];
    planetaryForce(theDomain->thePlanets+0, xyz, Fxyz);
    double Fp1 = cosp * Fxyz[1] - sinp * Fxyz[0];

    planetaryForce(theDomain->thePlanets+1, xyz, Fxyz);
    double Fp2 = cosp * Fxyz[1] - sinp * Fxyz[0];

    if(r > 1.0)
    {
        Q[0] = prim[RHO] * r * Fp1;
        Q[1] = prim[RHO] * r * Fp2;
        Q[2] = 0.0;
        Q[3] = 0.0;
    }
    else
    {
        Q[0] = 0.0;
        Q[1] = 0.0;
        Q[2] = prim[RHO] * r * Fp1;
        Q[3] = prim[RHO] * r * Fp2;
    }

    Q[4] = rho;
    Q[5] = rho * r*cosp;
    Q[6] = rho * r*sinp;
    Q[7] = rho * r*r*cos2p;
    Q[8] = rho * r*r*sin2p;

}

