#ifndef VDISCO_VUTIL_H
#define VDUSCO_VUTIL_H

#include <hdf5.h>

void getH5dims(char *file, char *group, char *dset, hsize_t *dims);
void readSimple(char *file, char *group, char *dset, void *data, hid_t type);
void readString(char *file, char *group, char *dset, char *buf, int len);
void readPatch(char *file, char *group, char *dset, void *data, hid_t type,
               int dim, int *start, int *loc_size, int *glo_size);

#endif
