#include "../paul.h"
#include "../geometry.h"
#include "../planet.h"


/*
 *  No additional reports at all!
 */

void setReportParams(struct domain *theDomain)
{
    // Silence is golden.
}


int num_shared_reports()
{
    return 0;
}

int num_distributed_aux_reports()
{
    return 0;
}

int num_distributed_integral_reports()
{
    return 0;
}

void get_shared_reports(double *Q, struct domain *theDomain)
{
    // Silence is golden.
}

void get_distributed_aux_reports(double *Q, struct domain *theDomain)
{
    // Silence is golden.
}

void get_distributed_integral_reports(const double *x, const double *prim,
                                      double *Q, struct domain *theDomain)
{
    // Silence is golden.
}

