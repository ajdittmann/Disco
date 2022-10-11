#ifndef DISCO_ANALYSIS_H
#define DISCO_ANALYSIS_H

#include "paul.h"

void setDiagParams(struct domain *theDomain);
int num_diagnostics(void);
int num_snapshot_rz(void);
int num_snapshot_arr(void);

void get_diagnostics(const double *x, const double *prim, double *Qrz, 
                        struct domain * theDomain );
void get_snapshot_rz(const double *x, const double *prim, double *Qrz, 
                        struct domain * theDomain );
void get_snapshot_arr( double * x , double * prim , double * Qarr, 
                           struct domain * theDomain );


void setup_diagnostics(struct domain * theDomain);
void setup_snapshot(struct domain * theDomain);

void zero_diagnostics( struct domain * theDomain );
void avg_diagnostics( struct domain * theDomain );
void adjust_RK_diag(struct domain *theDomain, double RK);
void copy_RK_diag(struct domain *theDomain);
void add_diagnostics( struct domain * theDomain , double dt );

void snapshot(struct domain *theDomain, char filename[]);

void free_diagnostics(struct diagnostic_avg *theTools);
void free_snapshot(struct snapshot_data *theSnap);

void zero_snapshot( struct domain * theDomain );
void calc_snapshot_rz(struct domain *theDomain);
void calc_snapshot_arr(struct domain *theDomain);

    
#endif
