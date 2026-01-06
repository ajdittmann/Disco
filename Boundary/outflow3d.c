
#include "../disco.h"
#include "../boundary.h"

void boundary_trans( struct domain * theDomain , int dim )
{
    if(dim == 1)
    {
        boundary_zerograd_rinn(theDomain, 1);
        boundary_zerograd_rout(theDomain, 1);
    }
    else if(dim == 2)
    {
        boundary_zerograd_zbot(theDomain, 1);
        boundary_zerograd_ztop(theDomain, 1);
    }
}

