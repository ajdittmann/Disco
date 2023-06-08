#include <stdio.h>
#include <math.h>
#include "calc.h"

// Local Functions

int bondi_newt(double Mdot, double GM, double gam, double rho0, 
                double *r, double *rho, double *u, double *P, int N)
{
    double qs = 0.25 * pow(2/(5-3*gam), (5-3*gam)/(2*gam-2));
    if(qs < 0.25 || fabs(3*gam-5) < 1.0e-6)
        qs = 0.25;
    else if (qs > 0.25*exp(1.5) || fabs(gam-1.0) < 1.0e-6)
        qs = 0.25*exp(1.5);
    double a0 = pow(4*M_PI*qs*GM*GM*rho0/Mdot, 1.0/3.0);
    double K = a0*a0 * pow(rho0, -gam+1.0) / gam;

    int i;
    for(i=0; i<N; i++)
    {
        double ui = bondi_newt_solve(Mdot, GM, r[i], gam, a0, K);

        rho[i] = -Mdot / (4*M_PI*r[i]*r[i] * ui);
        u[i] = ui;
        P[i] = K * pow(rho[i], gam);
    }

    return 0;
}

double bondi_newt_solve(double Mdot, double GM, double r, double gam, 
                        double a0, double K)
{
    double u0, u, umin, umax, u1, du, f, df;
    double TOL = 1.0e-6;
    int MAX_ITER = 100;

    double rhoFac = -Mdot / (4*M_PI*r*r);
    double a2Fac = gam * K;
    double Vr = -GM/r;
    double C = a0*a0/(gam-1.0);
    double rs = (5-3*gam)*GM/(4*a0*a0);
    double as = a0 * sqrt(2/(5-3*gam));

    if(r > rs)
    {
        umin = -as;
        umax = 0.0;
        u0 = -Mdot/(4*M_PI*r*r*pow((gam-1)*(C-Vr)/(K*gam), 1.0/(gam-1.0)));
    }
    else
    {
        umin = -1e100;
        umax = -as;
        u0 = -sqrt(2*(C-Vr));
    }

    //printf("%.6lg %.6lg %.6lg %.6lg\n", Mdot, GM, r, gam);
    //printf("%.6lg %.6lg %.6lg\n", a0, K);

    u1 = u0;

    int i = 0;
    do{
        u = u1;
        double rho = rhoFac / u;
        double drho = -rho / u;
        double a2 = a2Fac * pow(rho, gam-1.0);
        double da2 = (gam-1.0) * a2 / rho * drho;

        f = 0.5*u*u + a2/(gam-1.0) + Vr - C;
        df = u + da2/(gam-1.0);

        du = -f/df;

        u1 = u + du;

        if(u1 < umin)
            u1 = 0.5*(umin+u);
        else if(u1 > umax)
            u1 = 0.5*(umax+u);

        //printf("%d: u du = %.6lg %.6lg (%.6lg %.6lg)\n", i, u, du, f, df);

        i++;

    }while(fabs(du) > TOL && i < MAX_ITER);

    return u1;
}

int bondi_rel(double Mdot, double GM, double gam, double a0, 
                double *r, double *rho, double *u, double *P, int N)
{
    double as = bondi_rel_as(a0, gam);

    double us2 = as*as/(1+3*as*as);
    double rs = GM / (2*us2);
    double us = -sqrt(us2);
    double rhos = -Mdot / (4*M_PI*rs*rs*us);
    double K = as*as / ((1 - as*as/(gam-1)) * pow(rhos,gam-1) * gam);

    printf("M: %.8lg rs: %.8lg us: %.8lg as: %.8lg\n", GM, rs, us, as);

    double rho0 = pow(a0*a0 / (gam*K*(1.0-a0*a0/(gam-1))), 1.0/(gam-1));

    double h0 = 1.0 + gam/(gam-1) * K * pow(rho0, gam-1);

    printf("rho0: %.8lg h0: %.8lg k: %.8lg\n", rho0, h0, K);

    int i;
    for(i=0; i<N; i++)
    {
        double ui = bondi_rel_solve(Mdot, GM, r[i], gam, h0, K, as);

        rho[i] = -Mdot / (4*M_PI*r[i]*r[i] * ui);
        u[i] = ui;
        P[i] = K * pow(rho[i], gam);
    }

    return 0;
}

double bondi_rel_as(double a0, double gam)
{
    int MAX_ITER = 100;
    double TOL = 1.0e-13;
    double C = (1-a0*a0/(gam-1))*(1-a0*a0/(gam-1));
    double as2, as20, as21, das2;

    as20 = a0*a0;

    as21 = as20;

    double f, df;

    int i = 0;
    do
    {
        as2 = as21;
        f = (1+3*as2)*(1-as2/(gam-1))*(1-as2/(gam-1)) - C;
        df = (3*(1-as2/(gam-1)) - 2*(1+3*as2)/(gam-1))*(1-as2/(gam-1));

        das2 = -f/df;

        as21 = as2 + das2;

        //printf("%d: %.6lg %.6lg (%.6lg %.6lg)\n", i, as2, das2, f, df);

        i++;
    }while(fabs(das2) > TOL && i < MAX_ITER);

    double as = sqrt(as21);

    /*
    printf("N:%d as2=%.6lg as=%.6lg as*as=%.6lg aoo=%.6lg LHS: %.6lg RHS: %.6lg\n",
            i, as21, as, as*as, a0, 
            (1+3*as21)*(1-as21/(gam-1))*(1-as21/(gam-1)),
            (1-a0*a0/(gam-1))*(1-a0*a0/(gam-1)));
            */

    return as;
}

double bondi_rel_solve(double Mdot, double GM, double r, double gam, 
                        double h0, double K, double as)
{
    double TOL = 1.0e-12;
    int MAX_ITER = 100;

    double C = h0*h0;

    double us = -as/sqrt(1+3*as*as); 
    double rs = GM / (2*us*us);

    double u0, u1, u, du, f, df;
    double umin, umax;

    if(r > rs)
    {
        //u0 = -Mdot/(4*M_PI*r*r*pow((gam-1)*(C+2*GM/r)/(K*gam), 1.0/(gam-1.0)));
        double c1 = sqrt(C/(1-2*GM/r)) - 1.0;
        double a21 = (gam-1) * c1/(1+c1);
        double rho1 = pow(a21 / (K*gam*(1-a21/(gam-1))), 1.0/(gam-1));
        u0 = -Mdot/(4*M_PI*r*r*rho1);
        umin = us;
        umax = 0.0;
    }
    else
    {
        u0 = -sqrt(2*GM/r);
        umin = -1.0e100;
        umax = us;
    }
    u1 = u0;

    int i = 0;

    do
    {
        u = u1;

        double rho = -Mdot / (4*M_PI*r*r*u);
        double drho = -rho / u;
        double h = 1.0 + gam*K*pow(rho, gam-1)/(gam-1);
        double dh = gam*K*pow(rho, gam-2)*drho;

        f = h*h * (1 - 2*GM/r + u*u) - C;
        df = 2*h*(1 - 2*GM/r + u*u)*dh + 2*h*h*u;

        du = -f/df;

        u1 = u+du;

        if(u1 < umin)
            u1 = 0.5*(u+umin);
        if(u1 > umax)
            u1 = 0.5*(u+umax);
        
        //printf("%d: u du = %.6lg %.6lg (%.6lg %.6lg)\n", i, u, du, f, df);

        i++;
    }while(fabs(du) > TOL && i < MAX_ITER);

    //printf("%d: u du = %.6lg %.6lg (%.6lg %.6lg)\n", i, u, du, f, df);

    return u1;
}
