#include "../../disco.h"
#include "../metric.h"

/*
 * This is the Schwarzschild metric in Kerr-Schild spherical coordinates 
 * (t, r, phi, theta). These are related to the usual spherical Schwarzschild 
 * coordinates (t_sc, r, theta, phi) by: t = t_sc + 2M log|r/2M-1|, 
 */

static double om = 0.0; 
static double M = 0.0; 

void setMetricParams(struct domain *theDomain)
{
   om = theDomain->theParList.metricPar1;
   M = theDomain->theParList.metricPar2;
}

double metric_lapse(const double x[3])
{
    double r = x[0];
    return 1.0/sqrt(1.0 + 2*M/r);
}

void metric_shift(const double x[3], double b[3])
{
    double r = x[0];
    double a2 = 1.0 / (1.0 + 2*M/r);

    b[0] = 2*M/r * a2;
    b[1] = 0.0;
    b[2] = 0.0;
}

void metric_gam(const double x[3], double gam[9])
{
    double r = x[0];
    double th = x[2];


    double H = 2*M/r;
    double st = sin(th);

    gam[0] = 1.0 + H;
    gam[1] = 0.0;
    gam[2] = 0.0;
    gam[3] = 0.0;
    gam[4] = r*r*st*st;
    gam[5] = 0.0;
    gam[6] = 0.0;
    gam[7] = 0.0;
    gam[8] = r*r;
}

void metric_igam(const double x[3], double igam[9])
{
    double r = x[0];
    double th = x[2];

    double H = 2*M/r;
    double st = sin(th);

    igam[0] = 1.0 / (1.0+H);
    igam[1] = 0.0;
    igam[2] = 0.0;
    igam[3] = 0.0;
    igam[4] = 1.0/(r*r*st*st);
    igam[5] = 0.0;
    igam[6] = 0.0;
    igam[7] = 0.0;
    igam[8] = 1.0/(r*r);
}

double metric_jacobian(const double x[3])
{
    return x[0]*x[0]*sin(x[2]);
}

void metric_der_g(const double x[3], int i, double dg[16])
{
    int mu;
    for(mu=0; mu<16; mu++)
        dg[mu] = 0.0;

    if(i != 1 && i != 3)
        return;

    double r = x[0];
    double th = x[2];

    double st = sin(th);
    double ct = cos(th);

    double MoR2 = M/(r*r);
    double MoR3 = M/(r*r*r);
    double MoR4 = M/(r*r*r*r);

    if(i == 1)
    {
        dg[0]  = -2*MoR2;       // 00
        dg[1]  = -2*MoR2;       // 01
        dg[5]  =  -2*MoR2;      // 11
        dg[10] =  2*r*st*st;    // 22
        dg[15] = 2*r;           // 33
        dg[4]  = dg[1];
    }
    else if(i == 3)
    {
        dg[10] = 2*r*r*st*ct;   // 33
    }
}

void metric_der_lapse(const double x[3], double da[4])
{
    double r = x[0];

    double MoR2 = M/(r*r);
    double a = 1.0/sqrt(1.0+2*M/r);

    da[0] = 0.0;
    da[1] = a*a*a * MoR2;
    da[2] = 0.0;
    da[3] = 0.0;
}

void metric_der_shift(const double x[3], double db[12])
{
    double r = x[0];

    int i;

    for(i=0; i<12; i++)
        db[i] = 0.0;

    //dr b^r
    db[3] = - 2*M / ((r+2*M)*(r+2*M));
}

int metric_killing(int mu)
{
    if(mu == 1 || mu == 3)
        return 0;
    return 1;
}
