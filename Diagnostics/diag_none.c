
#include "../paul.h"
#include "../geometry.h"
#include "../planet.h"


void setDiagParams( struct domain * theDomain )
{
    // Silence is golden.
}

int num_diagnostics(void)
{
    return 0;
}

int num_snapshot_rz(void)
{
    return 0;
}

int num_snapshot_arr(void)
{
    return 0;
}


void get_diagnostics(const double *x, const double *prim, double *Qrz, 
                        struct domain * theDomain )
{
    // Silence is golden.
}


void get_snapshot_rz(const double *x, const double *prim, double *Qrz, 
                        struct domain * theDomain )
{
    // Silence is golden.
}

void get_snapshot_arr(const double *x, const double *prim, double *Qarr, 
                        struct domain * theDomain )
{
    // Silence is golden.
}
