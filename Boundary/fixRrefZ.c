
#include "../disco.h"
#include "../boundary.h"

void boundary_trans( struct domain * theDomain , int dim )
{
    if(dim == 1)
    {
        boundary_fixed_rinn(theDomain);
        boundary_fixed_rout(theDomain);
    }
    else if(dim == 2)
    {
        boundary_reflect_zbot(theDomain);
        boundary_reflect_ztop(theDomain);
    }
}

