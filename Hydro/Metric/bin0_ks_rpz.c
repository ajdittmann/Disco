#include "../../disco.h"
#include "../metric.h"

/*
 * This is an approximate metric for binary black holes, constructed as the
 * superposition of two kerr-schild metrics.  It is the initial condition
 * for some iterative schemes that fully solve the constraint equations to
 * generate initial data for numerical relativity codes.
 * See: Marronnetti et al. 2000 (gr-qc/0001077v2)
 *      Yo et al. 2004 (PRD.10.084033)
 *
 * An improvement would be to incorporate:
 * Tichy et al. 2002 (gr-qc/0207011) and/or
 * Nissanke 2006 (gr-qc/0509128)
 *
 * And finally to use:
 * Ireland et al. 2016 (gr-qc/1512.05650)
 */

static double M = 0.0; 
static double q = 1.0;
static double a = 10.0;
static double M1;
static double M2;
static double a1;
static double a2;
static double om;

void setMetricParams(struct domain *theDomain)
{
   M = theDomain->theParList.metricPar2;
   
    M1 = M/(1+q);
    M2 = q*M/(1+q);
    a1 = q*a/(1+q);
    a2 = -a/(1+q);
    om = sqrt(M/(a*a*a));
}

double metric_lapse(const double X[3])
{
    double x = X[0]*cos(X[1]);
    double y = X[0]*sin(X[1]);
    double z = X[2];

    double R1 = sqrt((x-a1)*(x-a1) + y*y + z*z);
    double R2 = sqrt((x-a2)*(x-a2) + y*y + z*z);

    double al1 = 1.0/sqrt(1.0 + 2*M1/R1);
    double al2 = 1.0/sqrt(1.0 + 2*M2/R2);

    return al1*al2;
}

void metric_shift(const double X[3], double b[3])
{
    double r = X[0];
    double cosp = cos(X[1]);
    double sinp = sin(X[1]);

    double x = r * cosp;
    double y = r * sinp;
    double z = X[2];

    double R1 = sqrt((x-a1)*(x-a1) + y*y + z*z);
    double R2 = sqrt((x-a2)*(x-a2) + y*y + z*z);

    double H1 = M1/R1;
    double H2 = M2/R2;

    double k1x = (a1-x)/R1;
    double k1y = -y/R1;
    double k1z = -z/R1;
    double k2x = (a2-x)/R2;
    double k2y = -y/R2;
    double k2z = -z/R2;

    double be1x = -2*H1*k1x/(1+2*H1);
    double be1y = -2*H1*k1y/(1+2*H1);
    double be1z = -2*H1*k1z/(1+2*H1);
    double be2x = -2*H2*k2x/(1+2*H2);
    double be2y = -2*H2*k2y/(1+2*H2);
    double be2z = -2*H2*k2z/(1+2*H2);
    

    b[0] = cosp*(be1x+be2x) + sinp*(be1y+be2y);
    b[1] = (-sinp*(be1x+be2x) + cosp*(be1y+be2y))/r + om;
    b[2] = be1z + be2z;
}

void metric_gam(const double X[3], double gam[9])
{
    double r = X[0];
    double cp = cos(X[1]);
    double sp = sin(X[1]);

    double x = r * cp;
    double y = r * sp;
    double z = X[2];

    double R1 = sqrt((x-a1)*(x-a1) + y*y + z*z);
    double R2 = sqrt((x-a2)*(x-a2) + y*y + z*z);

    double H1 = M1/R1;
    double H2 = M2/R2;

    //Covariant (w.r.t. Minkowski) k components in r/phi
    double k1r = (cp*a1 - r) / R1;
    double k1p = r*(-sp*a1) / R1;
    double k1z = -z/R1;
    double k2r = (cp*a2 - r) / R2;
    double k2p = r*(-sp*a2) / R2;
    double k2z = -z/R2;

    gam[0] = 1.0 + 2*H1*k1r*k1r + 2*H2*k2r*k2r;
    gam[1] = 2*H1*k1r*k1p + 2*H2*k2r*k2p;
    gam[2] = 2*H1*k1r*k1z + 2*H2*k2r*k2z;
    gam[3] = gam[1];
    gam[4] = r*r + 2*H1*k1p*k1p + 2*H2*k2p*k2p;
    gam[5] = 2*H1*k1p*k1z + 2*H2*k2p*k2z;
    gam[6] = gam[2];
    gam[7] = gam[5];
    gam[8] = 1.0 + 2*H1*k1z*k1z + 2*H2*k2z*k2z;
}

void metric_igam(const double X[3], double igam[9])
{
    double r = X[0];
    double cosp = cos(X[1]);
    double sinp = sin(X[1]);

    double x = r * cosp;
    double y = r * sinp;
    double z = X[2];

    double R1 = sqrt((x-a1)*(x-a1) + y*y + z*z);
    double R2 = sqrt((x-a2)*(x-a2) + y*y + z*z);

    double H1 = M1/R1;
    double H2 = M2/R2;

    double iR12 = 1.0/(R1*R1);
    double iR22 = 1.0/(R2*R2);

    double Px2 = y*y+z*z;

    double idetgam = 1.0/(1.0 + 2.0*H1 + 2.0*H2
                            + 4.0*H1*H2*a*a*Px2*iR12*iR22);

    double igxx = idetgam*(1.0 + 2*H1*Px2*iR12 + 2*H2*Px2*iR22);
    double igxy = idetgam*(2*H1*(a1-x)*y*iR12 + 2*H2*(a2-x)*y*iR22);
    double igxz = idetgam*(2*H1*(a1-x)*z*iR12 + 2*H2*(a2-x)*z*iR22);
    double igyy = idetgam*(1 + 2*H1*(R1*R1-y*y)*iR12 + 2*H2*(R2*R2-y*y)*iR22
                            + 4*H1*H2*a*a*z*z*iR12*iR22);
    double igyz = idetgam*(-2*H1*y*z*iR12 - 2*H2*y*z*iR22
                            - 4*H1*H2*a*a*y*z*iR12*iR22);
    double igzz = idetgam*(1 + 2*H1*(R1*R1-z*z)*iR12 + 2*H2*(R2*R2-z*z)*iR22
                            + 4*H1*H2*a*a*y*y*iR12*iR22);

    
    igam[0] = cosp*cosp*igxx + 2*cosp*sinp*igxy + sinp*sinp*igyy;
    igam[1] = (cosp*sinp*(-igxx+igyy) + (cosp*cosp-sinp*sinp)*igxy) / r;
    igam[2] = cosp*igxz + sinp*igyz;
    igam[3] = igam[1];
    igam[4] = (sinp*sinp*igxx - 2*cosp*sinp*igxy + cosp*cosp*igyy) / (r*r);
    igam[5] = (-sinp*igxz + cosp*igyz)/r;
    igam[6] = igam[2];
    igam[7] = igam[5];
    igam[8] = igzz;
}

double metric_jacobian(const double X[3])
{
    double r = X[0];
    double cosp = cos(X[1]);
    double sinp = sin(X[1]);

    double x = r * cosp;
    double y = r * sinp;
    double z = X[2];

    double R1 = sqrt((x-a1)*(x-a1) + y*y + z*z);
    double R2 = sqrt((x-a2)*(x-a2) + y*y + z*z);

    double H1 = M1/R1;
    double H2 = M2/R2;

    double iR12 = 1.0/(R1*R1);
    double iR22 = 1.0/(R2*R2);

    double Px2 = y*y+z*z;

    double detgam = 1.0 + 2.0*H1 + 2.0*H2 + 4.0*H1*H2*a*a*Px2*iR12*iR22;
    double al1 = 1.0/sqrt(1.0 + 2*M1/R1);
    double al2 = 1.0/sqrt(1.0 + 2*M2/R2);
    double al = al1*al2;

    return al * detgam * r;
}

void metric_der_g(const double X[3], int i, double dg[16])
{
    int mu;
    for(mu=0; mu<16; mu++)
        dg[mu] = 0.0;

    if(i == 0)
        return;

    double r = X[0];
    double cp = cos(X[1]);
    double sp = sin(X[1]);

    double x = r * cp;
    double y = r * sp;
    double z = X[2];

    double R1 = sqrt((x-a1)*(x-a1) + y*y + z*z);
    double R2 = sqrt((x-a2)*(x-a2) + y*y + z*z);

    double H1 = M1/R1;
    double H2 = M2/R2;

    double iR12 = 1.0/(R1*R1);
    double iR22 = 1.0/(R2*R2);

    //Covariant (w.r.t. Minkowski) k components in r/phi
    double k1r = (cp*a1 - r) / R1;
    double k1p = r*(-sp*a1) / R1;
    double k1z = -z/R1;
    double k2r = (cp*a2 - r) / R2;
    double k2p = r*(-sp*a2) / R2;
    double k2z = -z/R2;

    double al1 = 1.0/sqrt(1.0 + 2*M1/R1);
    double al2 = 1.0/sqrt(1.0 + 2*M2/R2);

    double dal1, dal2, dr2, dH1, dH2, dk1r, dk1p, dk1z, dk2r, dk2p, dk2z;

    if(i == 1) // d/dr
    {
        dr2 = 2*r;
        double dR1 = (r - cp*a1) / R1;
        double dR2 = (r - cp*a2) / R2;

        dal1 = -0.5 * al1*al1*al1 * (-2*M1*iR12*dR1);
        dal2 = -0.5 * al2*al2*al2 * (-2*M2*iR22*dR2);

        dH1 = -M1 * iR12 * dR1;
        dH2 = -M2 * iR22 * dR2;

        dk1r = (-R1-(cp*a1-r)*dR1) * iR12;
        dk2r = (-R2-(cp*a2-r)*dR2) * iR22;
        dk1p = -sp*a1 * (R1 - r*dR1) * iR12;
        dk2p = -sp*a2 * (R2 - r*dR2) * iR22;
        dk1z = z * iR12 * dR1;
        dk2z = z * iR22 * dR2;
    }
    else if(i == 2) // d/dphi
    {
        dr2 = 0.0;
        double dR1 = r*a1*sp / R1;
        double dR2 = r*a2*sp / R2;

        dal1 = -0.5 * al1*al1*al1 * (-2*M1*iR12*dR1);
        dal2 = -0.5 * al2*al2*al2 * (-2*M2*iR22*dR2);

        dH1 = -M1 * iR12 * dR1;
        dH2 = -M2 * iR22 * dR2;

        dk1r = (-sp*a1*R1-(cp*a1-r)*dR1) * iR12;
        dk2r = (-sp*a2*R2-(cp*a2-r)*dR2) * iR22;
        dk1p = -r*a1 * (cp*R1 - sp*dR1) * iR12;
        dk2p = -r*a2 * (cp*R2 - sp*dR2) * iR22;
        dk1z = z * iR12 * dR1;
        dk2z = z * iR22 * dR2;
    }
    else if(i == 3) // d/dz
    {
        dr2 = 0.0;
        double dR1 = z / R1;
        double dR2 = z / R2;

        dal1 = -0.5 * al1*al1*al1 * (-2*M1*iR12*dR1);
        dal2 = -0.5 * al2*al2*al2 * (-2*M2*iR22*dR2);

        dH1 = -M1 * iR12 * dR1;
        dH2 = -M2 * iR22 * dR2;

        dk1r = -(cp*a1-r) * iR12 * dR1;
        dk2r = -(cp*a2-r) * iR22 * dR2;
        dk1p = r*sp*a1 * iR12 * dR1;
        dk2p = r*sp*a2 * iR22 * dR2;
        dk1z = -(R1 - z*dR1) * iR12;
        dk2z = -(R2 - z*dR2) * iR22;
    }

    double al = al1*al2;
    double dal = al1*dal2 + dal1*al2;

    double be1r = -2*H1*k1r/(1+2*H1);
    double be1p = -2*H1*k1p/(r*r*(1+2*H1));
    double be1z = -2*H1*k1z/(1+2*H1);
    double be2r = -2*H2*k2r/(1+2*H2);
    double be2p = -2*H2*k2p/(r*r*(1+2*H2));
    double be2z = -2*H2*k2z/(1+2*H2);

    double dbe1r = -2*(dH1*k1r + H1*dk1r*(1+2*H1)) / ((1+2*H1)*(1+2*H1));
    double dbe1p = -2*(dH1*k1p + H1*dk1p*(1+2*H1)) / (r*r*(1+2*H1)*(1+2*H1)) - be1p*dr2/(r*r);
    double dbe1z = -2*(dH1*k1z + H1*dk1z*(1+2*H1)) / ((1+2*H1)*(1+2*H1));
    double dbe2r = -2*(dH2*k2r + H2*dk2r*(1+2*H2)) / ((1+2*H2)*(1+2*H2));
    double dbe2p = -2*(dH2*k2p + H2*dk2p*(1+2*H2)) / (r*r*(1+2*H2)*(1+2*H2)) - be2p*dr2/(r*r);
    double dbe2z = -2*(dH2*k2z + H2*dk2z*(1+2*H2)) / ((1+2*H2)*(1+2*H2));

    double ber = be1r+be2r;
    double bep = be1p+be2p+om;
    double bez = be1z+be2z;
    double dber = dbe1r+dbe2r;
    double dbep = dbe1p+dbe2p;
    double dbez = dbe1z+dbe2z;

    double grr = 1.0 + 2*H1*k1r*k1r + 2*H2*k2r*k2r;
    double grp = 2*H1*k1r*k1p + 2*H2*k2r*k2p;
    double grz = 2*H1*k1r*k1z + 2*H2*k2r*k2z;
    double gpp = r*r + 2*H1*k1p*k1p + 2*H2*k2p*k2p;
    double gpz = 2*H1*k1p*k1z + 2*H2*k2p*k2z;
    double gzz = 1.0 + 2*H1*k1z*k1z + 2*H2*k2z*k2z;

    //Spatially covariant shift beta_i = g_0i
    double belr = grr*ber + grp*bep + grz*bez;
    double belp = grp*ber + gpp*bep + gpz*bez;
    double belz = grz*ber + gpz*bep + gzz*bez;
    
    double dgrr, dgrp, dgrz, dgpp, dgpz, dgzz;

    dgrr = 2*dH1*k1r*k1r + 4*H1*k1r*dk1r + 2*dH2*k2r*k2r + 4*H2*k2r*dk2r;
    dgrp = 2*dH1*k1r*k1p + 2*H1*dk1r*k1p + 2*H1*k1r*dk1p
             + 2*dH2*k2r*k2p + 2*H2*dk2r*k2p + 2*H2*k2r*dk2p;
    dgrz = 2*dH1*k1r*k1z + 2*H1*dk1r*k1z + 2*H1*k1r*dk1z
             + 2*dH2*k2r*k2z + 2*H2*dk2r*k2z + 2*H2*k2r*dk2z;
    dgpp = dr2 + 2*dH1*k1p*k1p + 4*H1*k1p*dk1p + 2*dH2*k2p*k2p
                + 4*H2*k2p*dk2p;
    dgpz = 2*dH1*k1p*k1z + 2*H1*dk1p*k1z + 2*H1*k1p*dk1z
                + 2*dH2*k2p*k2z + 2*H2*dk2p*k2z + 2*H2*k2p*dk2z;
    dgzz = 2*dH1*k1z*k1z + 4*H1*k1z*dk1z + 2*dH2*k2z*k2z + 4*H2*k2z*dk2z;

    double dbelr = dgrr*ber+grr*dber + dgrp*bep+grp*dbep + dgrz*bez+grz*dbez;
    double dbelp = dgrp*ber+grp*dber + dgpp*bep+gpp*dbep + dgpz*bez+gpz*dbez;
    double dbelz = dgrz*ber+grz*dber + dgpz*bep+gpz*dbep + dgzz*bez+gzz*dbez;

    double db2 = dber*belr+ber*dbelr + dbep*belp+bep*dbelp
                + dbez*belz+bez*dbelz;

    dg[0] = -2*al*dal + db2;
    dg[1] = dbelr;
    dg[2] = dbelp;
    dg[3] = dbelz;
    dg[4] = dg[1];
    dg[5] = dgrr;
    dg[6] = dgrp;
    dg[7] = dgrz;
    dg[8] = dg[2];
    dg[9] = dgrp;
    dg[10] = dgpp;
    dg[11] = dgpz;
    dg[12] = dg[3];
    dg[13] = dgrz;
    dg[14] = dgpz;
    dg[15] = dgzz;
}

void metric_der_lapse(const double X[3], double da[4])
{
    double r = X[0];
    double cp = cos(X[1]);
    double sp = sin(X[1]);

    double x = r * cp;
    double y = r * sp;
    double z = X[2];

    double R1 = sqrt((x-a1)*(x-a1) + y*y + z*z);
    double R2 = sqrt((x-a2)*(x-a2) + y*y + z*z);

    double H1 = M1/R1;
    double H2 = M2/R2;

    double iR1 = 1.0/R1;
    double iR2 = 1.0/R2;

    double al1 = 1.0/sqrt(1.0 + 2*H1);
    double al2 = 1.0/sqrt(1.0 + 2*H2);

    double dR1dr, dR1dp, dR1dz, dR2dr, dR2dp, dR2dz;

    dR1dr = (r - cp*a1) * iR1;
    dR2dr = (r - cp*a2) * iR2;
    dR1dp = r*a1*sp * iR1;
    dR2dp = r*a2*sp * iR2;
    dR1dz = z * iR1;
    dR2dz = z * iR2;

    double al13 = al1*al1*al1;
    double al23 = al2*al2*al2;
    double M1_o_R12 = M1 * iR1*iR1;
    double M2_o_R22 = M2 * iR2*iR2;

    da[0] = 0.0;
    da[1] = al1 * al23*M2_o_R22*dR2dr + al13*M1_o_R12*dR1dr * al2;
    da[2] = al1 * al23*M2_o_R22*dR2dp + al13*M1_o_R12*dR1dp * al2;
    da[3] = al1 * al23*M2_o_R22*dR2dz + al13*M1_o_R12*dR1dz * al2;
}

void metric_der_shift(const double X[3], double db[12])
{
    double r = X[0];
    double cp = cos(X[1]);
    double sp = sin(X[1]);

    double x = r * cp;
    double y = r * sp;
    double z = X[2];

    double R1 = sqrt((x-a1)*(x-a1) + y*y + z*z);
    double R2 = sqrt((x-a2)*(x-a2) + y*y + z*z);

    double H1 = M1/R1;
    double H2 = M2/R2;

    double iR1 = 1.0/R1;
    double iR2 = 1.0/R2;
    double M1_o_R12 = M1 / (R1*R1);
    double M2_o_R22 = M2 / (R2*R2);
    
    //Contravariant (w.r.t. Minkowski) k components in r/phi
    double k1[3], k2[3];
    k1[0] = (cp*a1 - r) * iR1;
    k1[1] = (-sp*a1) * iR1 / r;
    k1[2] = -z * iR1;
    k2[0] = (cp*a2 - r) * iR2;
    k2[1] = (-sp*a2) * iR2 / r;
    k2[2] = -z * iR2;

    double dR1[3], dR2[3], dH1[3], dH2[3];

    dR1[0] = (r - cp*a1) * iR1;
    dR1[1] = r*a1*sp * iR1;
    dR1[2] = z * iR1;
    dR2[0] = (r - cp*a2) * iR2;
    dR2[1] = r*a2*sp * iR2;
    dR2[2] = z * iR2;

    int i;
    for(i=0; i<3; i++)
    {
        dH1[i] = -M1_o_R12 * dR1[i];
        dH2[i] = -M2_o_R22 * dR2[i];
    }

    double dk1[9], dk2[9];  //dk[3*i+j] = d_i k^j

    dk1[3*0+0] = (-R1 - (cp*a1-r)*dR1[0]) * iR1*iR1;
    dk1[3*0+1] = sp*a1 * (iR1*iR1*dR1[0]/r + iR1/(r*r));
    dk1[3*0+2] = z * iR1*iR1*dR1[0];
    dk1[3*1+0] = (-sp*a1*R1 - (cp*a1-r)*dR1[1]) * iR1*iR1;
    dk1[3*1+1] = -a1*(cp*R1 - sp*dR1[1]) * iR1*iR1 / r;
    dk1[3*1+2] = z * iR1*iR1*dR1[1];
    dk1[3*2+0] = -(cp*a1-r)*iR1*iR1*dR1[2];
    dk1[3*2+1] = sp*a1 * iR1*iR1 * dR1[2] / r;
    dk1[3*2+2] = -(R1 - z*dR1[2]) * iR1*iR1;

    dk2[3*0+0] = (-R2 - (cp*a2-r)*dR2[0]) * iR2*iR2;
    dk2[3*0+1] = sp*a2 * (iR2*iR2*dR2[0]/r + iR2/(r*r));
    dk2[3*0+2] = z * iR2*iR2*dR2[0];
    dk2[3*1+0] = (-sp*a2*R2 - (cp*a2-r)*dR2[1]) * iR2*iR2;
    dk2[3*1+1] = -a2*(cp*R2 - sp*dR2[1]) * iR2*iR2 / r;
    dk2[3*1+2] = z * iR2*iR2*dR2[1];
    dk2[3*2+0] = -(cp*a2-r)*iR2*iR2*dR2[2];
    dk2[3*2+1] = sp*a2 * iR2*iR2 * dR2[2] / r;
    dk2[3*2+2] = -(R2 - z*dR2[2]) * iR2*iR2;

    double denom1 = 1.0 / ((1+2*H1)*(1+2*H1));
    double denom2 = 1.0 / ((1+2*H2)*(1+2*H2));

    int j;
    for(j=0; j<3; j++)
        db[j] = 0.0;

    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
        {
            db[3*(i+1)+j] = -2*(dH1[i]*k1[j] + H1*dk1[3*i+j]*(1+2*H1)) * denom1
                            - 2*(dH2[i]*k2[j] + H2*dk2[3*i+j]*(1+2*H2)) * denom2;
        }
}

int metric_killing(int mu)
{
    if(mu == 0)
        return 1;
    return 0;
}
