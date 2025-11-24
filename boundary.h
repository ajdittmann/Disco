#ifndef DISCO_BOUNDARY_H
#define DISCO_BOUNDARY_H

#include "disco.h"

// For initialization
void setBCParams(struct domain *theDomain);

//The boundary function, defined in Boundary/
void boundary_trans( struct domain * theDomain , int dim );

//Boundary options
void boundary_fixed_rinn( struct domain *theDomain);
void boundary_fixed_rout( struct domain *theDomain);
void boundary_fixed_zbot( struct domain *theDomain);
void boundary_fixed_ztop( struct domain *theDomain);

void boundary_zerograd_rinn( struct domain *theDomain, int diode);
void boundary_zerograd_rout( struct domain *theDomain, int diode);
void boundary_zerograd_zbot( struct domain *theDomain, int diode);
void boundary_zerograd_ztop( struct domain *theDomain, int diode);

void boundary_reflect_rinn( struct domain *theDomain);
void boundary_reflect_rout( struct domain *theDomain);
void boundary_reflect_zbot( struct domain *theDomain);
void boundary_reflect_ztop( struct domain *theDomain);

void boundary_fixed_q_rinn( struct domain *theDomain, int *q, int nq);
void boundary_fixed_q_rout( struct domain *theDomain, int *q, int nq);
void boundary_fixed_q_zbot( struct domain *theDomain, int *q, int nq);
void boundary_fixed_q_ztop( struct domain *theDomain, int *q, int nq);

void boundary_fixed_horizon( struct domain *theDomain);

void boundary_fixed_phi_rinn( struct domain *theDomain, double phia,
                                double phib);
void boundary_fixed_phi_rout( struct domain *theDomain, double phia,
                                double phib);
void boundary_fixed_phi_zbot( struct domain *theDomain, double phia,
                                double phib);
void boundary_fixed_phi_ztop( struct domain *theDomain, double phia,
                                double phib);

void boundary_noslip_rinn( struct domain *theDomain);
void boundary_noslip_rout( struct domain *theDomain);
void boundary_noslip_zbot( struct domain *theDomain);
void boundary_noslip_ztop( struct domain *theDomain);

//Low-level routines for boundary
void set_cell_init(struct cell *c, const double *r_jph, const double *z_kph,
                   int j, int k);
void set_cell_init_q(struct cell *c, const double *r_jph, const double *z_kph, 
                     int j, int k, int *qarr, int nq);
void set_cells_copy(struct cell *c, int Np, const struct face *theFaces, 
                    int n0, int n1, int LR);
void set_cells_copy_distant(struct cell *c, int Np, struct cell *c1, int Np1);

#endif
