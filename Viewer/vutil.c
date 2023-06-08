#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include "vutil.h"

void getH5dims(char *file, char *group, char *dset, hsize_t *dims)
{
    hid_t h5fil = H5Fopen( file , H5F_ACC_RDONLY , H5P_DEFAULT );
    hid_t h5grp = H5Gopen1( h5fil , group );
    hid_t h5dst = H5Dopen1( h5grp , dset );
    hid_t h5spc = H5Dget_space( h5dst );

    H5Sget_simple_extent_dims( h5spc , dims , NULL);

    H5Sclose( h5spc );
    H5Dclose( h5dst );
    H5Gclose( h5grp );
    H5Fclose( h5fil );
}

void readSimple(char *file, char *group, char *dset, void *data, hid_t type)
{
    hid_t h5fil = H5Fopen( file , H5F_ACC_RDONLY , H5P_DEFAULT );
    hid_t h5grp = H5Gopen1( h5fil , group );
    hid_t h5dst = H5Dopen1( h5grp , dset );

    H5Dread( h5dst , type , H5S_ALL , H5S_ALL , H5P_DEFAULT , data );

    H5Dclose( h5dst );
    H5Gclose( h5grp );
    H5Fclose( h5fil );
}

void readString(char *file, char *group, char *dset, char *buf, int len)
{
    hid_t strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, H5T_VARIABLE);

    hid_t h5fil = H5Fopen( file , H5F_ACC_RDONLY , H5P_DEFAULT );
    hid_t h5grp = H5Gopen1( h5fil , group );

    if(!H5Lexists(h5grp, dset, H5P_DEFAULT))
    {
        buf[0] = '\0';
        return;
    }
    hid_t h5dst = H5Dopen1( h5grp , dset );

    hid_t space = H5Dget_space(h5dst);
    hsize_t dims[1];
    H5Sget_simple_extent_dims(space, dims, NULL);

    char **rdata = (char **)malloc(dims[0] * sizeof(char *));
    H5Dread(h5dst, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata);
    strncpy(buf, rdata[0], len-1);
    buf[len-1] = '\0';

    H5Dvlen_reclaim(strtype, space, H5P_DEFAULT, rdata);
    free(rdata);
    H5Dclose(h5dst);
    H5Sclose(space);
    H5Tclose(strtype);
    H5Gclose(h5grp);
    H5Fclose(h5fil);
}

void readPatch(char *file, char *group, char *dset, void *data, hid_t type,
               int dim, int *start, int *loc_size, int *glo_size)
{
    hid_t h5fil = H5Fopen( file , H5F_ACC_RDONLY , H5P_DEFAULT );
    hid_t h5grp = H5Gopen1( h5fil , group );
    hid_t h5dst = H5Dopen1( h5grp , dset );

    hsize_t mdims[dim];
    hsize_t fdims[dim];

    hsize_t fstart[dim];
    hsize_t fstride[dim];
    hsize_t fcount[dim];
    hsize_t fblock[dim];

    int d;
    for( d=0 ; d<dim ; ++d )
    {
        mdims[d] = loc_size[d];
        fdims[d] = glo_size[d];

        fstart[d]  = start[d];
        fstride[d] = 1;
        fcount[d]  = loc_size[d];
        fblock[d]  = 1;
    }
    hid_t mspace = H5Screate_simple(dim,mdims,NULL);
    hid_t fspace = H5Screate_simple(dim,fdims,NULL);

    H5Sselect_hyperslab( fspace , H5S_SELECT_SET , fstart , fstride , fcount , fblock );

    H5Dread( h5dst , type , mspace , fspace , H5P_DEFAULT , data );

    H5Sclose( mspace );
    H5Sclose( fspace );
    H5Dclose( h5dst );
    H5Gclose( h5grp );
    H5Fclose( h5fil );
}
