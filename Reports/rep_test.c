#include "../disco.h"
#include "../geometry.h"
#include "../planet.h"

static double gamma_law = 0.0;
static int Npl = 0;

void setReportParams(struct domain *theDomain)
{
    gamma_law = theDomain->theParList.Adiabatic_Index;
    Npl = theDomain->Npl;
}

int num_shared_reports()
{
    return 0;
}

int num_distributed_aux_reports()
{
    return Npl * NUM_PL_AUX;
}

int num_distributed_integral_reports()
{
    return Npl;
}

void get_shared_reports(double *Q, struct domain *theDomain)
{
    //Silence is golden
}

void get_distributed_aux_reports(double *Q, struct domain *theDomain)
{
    int pq;

    for(pq=0; pq < Npl * NUM_PL_AUX; pq++)
        Q[pq] = theDomain->pl_aux[pq];
}

void get_distributed_integral_reports(double *x, double *prim, double *Q,
                                      struct domain *theDomain)
{
    int p;

    for(p=0; p<Npl; p++)
    {
        double xyz[3];
        get_xyz(x, xyz);
        Q[p] = prim[RHO] * planetaryPotential(theDomain->thePlanets+p, xyz);
    }
}

