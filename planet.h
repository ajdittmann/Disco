#ifndef DISCO_PLANET_H
#define DISCO_PLANET_H

#include "disco.h"

void setGravParams( struct domain * theDomain );

double phigrav( double M , double r , double eps , int type);
double fgrav( double M , double r , double eps , int type);
void adjust_gas( struct planet * pl , double * x , double * prim , double gam);
double planetaryPotential(struct planet *pl, const double *xyz);
void planetaryForce( struct planet * pl , const double *xyz, double *Fxyz);
void planet_src( struct planet * pl , const double * prim , double * cons ,
                const double * xp , const double * xm , const double *xyz,
                double dV, double dt, double *pl_gas_track);

void planet_init_kin(struct planet *pl, double *pl_kin);

void zeroAuxPlanets(struct domain *theDomain);
void setupPlanets(struct domain *theDomain);
void setPlanetsXYZ(struct domain *theDomain);
void initializePlanetTracking(struct domain *theDomain);
void updatePlanetsKinAux(struct domain *theDomain, double dt);
void movePlanetsLive(struct domain *theDomain);
void copyPlanetsRK( struct domain *theDomain);
void adjustPlanetsRKkin( struct domain *theDomain, double RK);
void adjustPlanetsRKaux( struct domain *theDomain, double RK);

// Functions in Planet/ setup
void setPlanetParams( struct domain * theDomain );
int planet_motion_analytic( void );
void initializePlanets( struct planet * thePlanets );
void movePlanets( struct planet * thePlanets , double t , double dt );
void forcePlanets( struct planet * thePlanets , double dt );

#endif
