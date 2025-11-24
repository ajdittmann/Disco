
#include "../disco.h"
#include "../boundary.h"

void boundary_trans( struct domain * theDomain , int dim )
{
    if(dim == 1)
    {
        boundary_reflect_rinn(theDomain);
        boundary_reflect_rout(theDomain);
    }
    else if(dim == 2)
    {
        boundary_reflect_zbot(theDomain);
        boundary_zerograd_ztop(theDomain, 1);
    }
}

