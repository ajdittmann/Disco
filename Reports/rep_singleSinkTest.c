#include "../paul.h"
#include "../geometry.h"
#include "../planet.h"

#define N_AUX_PER_PLANET 17

static double gamma_law = 0.0;
int Npl = 0;


/*
 * 
 */

void setReportParams(struct domain *theDomain)
{
    gamma_law = theDomain->theParList.Adiabatic_Index;
    Npl = theDomain->Npl;
}

int num_shared_reports()
{
    return 3;
}

int num_distributed_aux_reports()
{
    return N_AUX_PER_PLANET * Npl;
}

int num_distributed_integral_reports()
{
    return 0;
}

void get_shared_reports(double *Q, struct domain *theDomain)
{
    Q[0] = theDomain->thePlanets[0].M;
    Q[1] = theDomain->thePlanets[1].M;
    Q[2] = theDomain->thePlanets[2].M;
}

void get_distributed_aux_reports(double *Q, struct domain *theDomain)
{
    int idx_aux[N_AUX_PER_PLANET] = {
                PL_SNK_M, PL_GRV_JZ, PL_SNK_JZ, PL_SNK_SZ,
                PL_GRV_PX, PL_GRV_PY, PL_GRV_PZ,
                PL_SNK_PX, PL_SNK_PY, PL_SNK_PZ,
                PL_SNK_MX, PL_SNK_MY, PL_SNK_MZ,
                PL_GRV_K, PL_GRV_U, PL_SNK_K, PL_SNK_U};

    int p, q;
    for(p = 0; p < theDomain->Npl; p++)
    {
        int idxQ = p * N_AUX_PER_PLANET;
        int idxP = p * NUM_PL_AUX;
        for(q=0; q<N_AUX_PER_PLANET; q++)
            Q[idxQ + q] = theDomain->pl_aux[idxP + idx_aux[q]];
    }
}

void get_distributed_integral_reports(const double *x, const double *prim,
                                      double *Q, struct domain *theDomain)
{
    //Silence is golden.
}

