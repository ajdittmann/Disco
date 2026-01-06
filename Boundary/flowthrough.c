
#include "../disco.h"
#include "../boundary.h"


void boundary_trans( struct domain * theDomain , int dim )
{
    if(dim == 1)
    {
        boundary_noslip_rinn(theDomain);
        boundary_zerograd_rout(theDomain, 1);
        boundary_fixed_phi_rout(theDomain, 0.3*M_PI, 1.7*M_PI);
    }
    else if(dim == 2)
    {
        boundary_fixed_zbot(theDomain);
        boundary_fixed_ztop(theDomain);
    }
}

