#ifndef DISCO_HYDRO_H
#define DISCO_HYDRO_H

#include "disco.h"

void setHydroParams( struct domain * theDomain );
int set_B_flag();

double get_omega(const double *prim, const double *x);

void prim2cons(const double *prim, double *cons, const double *x,
               double dV, const double *xp, const double *xm);
void getUstar(const double *prim, double *Ustar, const double *x,
              double Sk, double Ss, const double *n, const double *Bpack);
void cons2prim(const double *cons, double *prim, const double *x, double dV,
               const double *xp, const double *xm);
void flux(const double *prim, double *flux, const double * x, const double * n,
          const double *xp, const double *xm);
void source(const double * prim , double * cons ,
            const double * xp , const double * xm , double dVdt);
void visc_flux(const double * prim, const double * gradr, const double * gradp,
               const double * gradz, double * flux,
               const double * x, const double * n);
void visc_source(const double * prim, const double * gradr, const double *gradp,
                 const double * gradz, double * cons, const double *xp,
                 const double *xm, double dVdt);
void prim_to_E(const double *prim, double *E, const double *x);
void flux_to_E( const double * Flux , const double * Ustr , const double * x, 
                double * E1_riemann , double * B1_riemann , 
                double * E2_riemann , double * B2_riemann , int dim );
void vel( const double * prim1 , const double * prim2 , 
         double * Sl , double * Sr , double * Ss , 
         const double * n , const double * x , double * Bpack );
double mindt(const double * prim , double w ,
             const double * xp , const double * xm );
void reflect_prims(double * prim, const double * x, int dim);
double bfield_scale_factor(double x, int dim);
#if ENABLE_CART_INTERP
double getCartInterpWeight(const double *x);
#endif

#endif
