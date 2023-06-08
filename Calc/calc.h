#ifndef DISCO_CALC
#define DISCO_CALC

int bondi_newt(double Mdot, double GM, double gam, double rho0, 
                double *r, double *rho, double *u, double *P, int N);
int bondi_rel(double Mdot, double GM, double gam, double a0, 
                double *r, double *rho, double *u, double *P, int N);

double bondi_newt_solve(double Mdot, double GM, double r, double gam, 
                        double a0, double K);
double bondi_rel_as(double a0, double gam);
double bondi_rel_solve(double Mdot, double GM, double r, double gam, 
                        double h0, double K, double as);

double magnetosonic_cf_int_newt(double rho, double rho0, double cs0,
                                    double cA0, double gam);
#endif
