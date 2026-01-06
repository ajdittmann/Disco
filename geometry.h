#ifndef DISCO_GEOMETRY_H
#define DISCO_GEOMETRY_H

#include "disco.h"

void setGeometryParams( struct domain * theDomain );

double get_dp( double phip , double phim );
double get_signed_dp( double phip , double phim );

double get_centroid(double xp, double xm, int dim);
double get_dL( const double * xp , const double * xm , int dim );
double get_dA( const double * xp , const double * xm , int dim );
double get_dV( const double * xp , const double * xm );
double get_scale_factor( const double * x, int dim);
double get_vol_element(const double *x);

void get_xyz(const double *x, double *xyz);
void get_rpz(const double *x, double *rpz);
void get_coords_from_xyz(const double *xyz, double *x);
void get_coords_from_rpz(const double *rpz, double *x);

void get_vec_rpz(const double *x, const double *v, double *vrpz);
void get_vec_from_rpz(const double *x, const double *vrpz, double *v);
void get_vec_xyz(const double *x, const double *v, double *vxyz);
void get_vec_from_xyz(const double *x, const double *vxyz, double *v);

void geom_grad(const double *prim, double *grad,
               const double *xp,const double *xm, double PLM, int dim, int LR);

#if ENABLE_CART_INTERP
void geom_interpolate(const double *prim, const double *gradr,
                      const double *gradp, const double *gradz,
                      const double *x, double dr, double dphi, double dz,
                      double * primI, double weight);

void geom_rebase_to_cart(const double *prim, const double *x,
                         double *cartPrim);
void geom_rebase_from_cart(const double *cartPrim, const double *x,
                           double *prim);
void geom_gradCart_to_grad(const double *cartGrad, const double *prim,
                           const double *x, double *grad, int dim);
void geom_cart_interp_grad_trans(const double *primL, const double *primR,
                                 const double *gradpL, const double *gradpR,
                                 double dpL, double dpR, const double *x,
                                 double dxL, double dxR, 
                                 double *gradCIL, double *gradCIR,
                                 int dim);
#endif

void get_centroid_arr(const double *xp, const double *xm, double *x);
void get_vec_contravariant(const double *x, double *v, double *vc);
void get_vec_covariant(const double *x, double *v, double *vc);

void geom_polar_vec_adjust(const double *xp, const double *xm, double *fac);

#endif
