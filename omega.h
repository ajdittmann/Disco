#ifndef DISCO_OMEGA_H
#define DISCO_OMEGA_H

#include "paul.h"

void setOmegaParams( struct domain * theDomain );

double mesh_om( const double *x);

double get_om( const double *x );
double get_om1( const double *x);
double get_om2( const double *x);
double get_height_om( const double *x);

double get_cs2( const double *x );
double get_nu( const double *x , const double *prim);

#endif
