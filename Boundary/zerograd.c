
#include "../disco.h"
#include "../boundary.h"

void boundary_trans( struct domain * theDomain , int dim )
{
    if(dim == 1)
    {
        boundary_zerograd_rinn(theDomain, 0);
        boundary_zerograd_rout(theDomain, 0);
    }
    else if(dim == 2)
    {
        boundary_zerograd_zbot(theDomain, 0);
        boundary_zerograd_ztop(theDomain, 0);
    }
}

