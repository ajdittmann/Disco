#include "paul.h"
#include "omega.h"
#include "geometry.h"

static int sinkRateType = 0;
static int sinkShapeType = 0;
static int sinkNumber = 0;
static double sinkPar1 = 0.0;
static double sinkPar2 = 0.0;
static double sinkPar3 = 0.0;
static double sinkPar4 = 0.0;
static double sinkPar5 = 0.0;
static int nozzleType;
static double nozzlePar1 = 0.0;
static double nozzlePar2 = 0.0;
static double nozzlePar3 = 0.0;
static double nozzlePar4 = 0.0;

static int coolType = 0;
static double coolPar1 = 0.0;
static double coolPar2 = 0.0;
static double coolPar3 = 0.0;
static double coolPar4 = 0.0;

static double gamma_law = 0.0;
static int twoD = 0;
static int Npl = 0;
static struct planet *thePlanets = NULL;

static double visc = 0.0;
static double Mach = 1.0;

static double rmax;
static double rmin;
static double zmax;
static double zmin;

static int DAMP_OUTER = 0;
static int DAMP_INNER = 0;
static int DAMP_LOWER = 0;
static int DAMP_UPPER = 0;
static double dampTimeInner = 0.0;
static double dampTimeOuter = 0.0;
static double dampTimeLower = 0.0;
static double dampTimeUpper = 0.0;
static double dampLenInner = 0.0;
static double dampLenOuter = 0.0;
static double dampLenUpper = 0.0;
static double dampLenLower = 0.0;

void planetaryForce( struct planet * , double , double , double , double * , double * , double * , int );
double phigrav( double , double , double , int);
void initial( double * , double * );
void prim2cons( double * , double * , double * , double );
double get_nu( const double *, const double *);

void setSinkParams(struct domain *theDomain)
{
    sinkShapeType = theDomain->theParList.sinkShapeType;
    sinkRateType = theDomain->theParList.sinkRateType;
    sinkNumber = theDomain->theParList.sinkNumber;
    sinkPar1 = theDomain->theParList.sinkPar1;
    sinkPar2 = theDomain->theParList.sinkPar2;
    sinkPar3 = theDomain->theParList.sinkPar3;
    sinkPar4 = theDomain->theParList.sinkPar4;
    sinkPar5 = theDomain->theParList.sinkPar5;
    nozzleType = theDomain->theParList.nozzleType;
    nozzlePar1 = theDomain->theParList.nozzlePar1;
    nozzlePar2 = theDomain->theParList.nozzlePar2;
    nozzlePar3 = theDomain->theParList.nozzlePar3;
    nozzlePar4 = theDomain->theParList.nozzlePar4;
    coolType = theDomain->theParList.coolType;
    coolPar1 = theDomain->theParList.coolPar1;
    coolPar2 = theDomain->theParList.coolPar2;
    coolPar3 = theDomain->theParList.coolPar3;
    coolPar4 = theDomain->theParList.coolPar4;

    rmax = theDomain->theParList.rmax;
    rmin = theDomain->theParList.rmin;
    zmax = theDomain->theParList.zmax;
    zmin = theDomain->theParList.zmin;

    gamma_law = theDomain->theParList.Adiabatic_Index;
    if(theDomain->Nz == 1)
        twoD = 1;

    visc = theDomain->theParList.viscosity;
    Mach = theDomain->theParList.Disk_Mach;

    thePlanets = theDomain->thePlanets;
    Npl = theDomain->Npl;

    DAMP_INNER = theDomain->theParList.dampInnerType;
    DAMP_OUTER = theDomain->theParList.dampOuterType;
    DAMP_UPPER = theDomain->theParList.dampUpperType;
    DAMP_LOWER = theDomain->theParList.dampLowerType;
    dampTimeInner = theDomain->theParList.dampTimeInner;
    dampLenInner = theDomain->theParList.dampLenInner;
    dampTimeOuter = theDomain->theParList.dampTimeOuter;
    dampLenOuter = theDomain->theParList.dampLenOuter;
    dampTimeUpper = theDomain->theParList.dampTimeUpper;
    dampLenUpper = theDomain->theParList.dampLenUpper;
    dampTimeLower = theDomain->theParList.dampTimeLower;
    dampLenLower = theDomain->theParList.dampLenLower;
}


void sink_src(double *prim, double *cons, double *xp, double *xm, double dV, double dt)
{
    if(nozzleType == 1)
    {
        double R0 = nozzlePar1;
        double dR = nozzlePar2;
        double v = nozzlePar3;
        double mach = nozzlePar4;

        double r = 0.5*(xp[0]+xm[0]);
        double phi = 0.5*(xp[1]+xm[1]);
        double z = 0.5*(xp[2]+xm[2]);
        if(phi > M_PI) phi -= 2*M_PI;

        if((r-R0)*(r-R0) + r*r*phi*phi + z*z > dR*dR)
            return;

        double Mdot = 1.0;
        double rhodot = Mdot * 3.0/(4.0*M_PI*dR*dR*dR);
        if(twoD)
            rhodot *= 2*sqrt(dR*dR - (r-R0)*(r-R0) - r*r*phi*phi);
        if(rhodot * dV*dt > cons[RHO])
            rhodot = cons[RHO] / dV*dt;

        double x[3] = {r, phi, z};
        double om = get_om(x);
        double cs2 = v*v / (mach*mach);
        double eps = cs2 / (gamma_law*(gamma_law-1));

        cons[RHO] += rhodot * dV*dt;
        cons[LLL] += rhodot * r*v * dV*dt;
        cons[TAU] += (0.5* rhodot * (v-r*om)*(v-r*om) + rhodot*eps) * dV*dt;
        if(NUM_Q > NUM_C)
            cons[NUM_C] += rhodot*dV*dt;
    }



    if(sinkShapeType != 0){

      double x[3];
      get_centroid_arr(xp, xm, x);
      double r = x[0];
      double phi = x[1];

      double rho = prim[RHO];
      double vr  = prim[URR];
      double vp  = prim[UPP]*r;
      double vz  = prim[UZZ];
      double press  = prim[PPP];
      double specenth = press*(1.0 + 1.0/(gamma_law - 1.0))/rho;

      double cosg = cos(phi);
      double sing = sin(phi);
      double gx = r*cosg;
      double gy = r*sing;

      double px, py, dx, dy, mag, eps, epsfactor, mag2;
      double surfdiff;
      int pi;
      int numSinks = Npl;

      double rbin = 1.0;
      double dxbin = 1.0;
      double dybin = 1.0;
      if (Npl == 2) {
        double cosp1 = cos(thePlanets[0].phi);
        double sinp1 = sin(thePlanets[0].phi);
        double rp1 = thePlanets[0].r;
        double px1 = rp1*cosp1;
        double py1 = rp1*sinp1;

        double cosp2 = cos(thePlanets[1].phi);
        double sinp2 = sin(thePlanets[1].phi);
        double rp2 = thePlanets[1].r;
        double px2 = rp2*cosp2;
        double py2 = rp2*sinp2;


        dxbin = px2-px1;
        dybin = py2-py1;
        rbin = sqrt( dxbin*dxbin + dybin*dybin );
      }
      if ((sinkNumber>0) && (sinkNumber<numSinks)) numSinks = sinkNumber;
      for (pi=0; pi<numSinks; pi++){
          double cosp = cos(thePlanets[pi].phi);
          double sinp = sin(thePlanets[pi].phi);
          px = thePlanets[pi].r*cosp;
          py = thePlanets[pi].r*sinp;

          dx = gx-px;
          dy = gy-py;
          mag2 = dx*dx + dy*dy;
          mag = sqrt(mag2);

          //the part that depends on sink profile/shape
          double arg = 0.0;
          if(sinkShapeType == 4){	//constant globally - use with care
            arg = 1.0;
          }
          else if(sinkShapeType == 3){	//constant if within sink radius
            if (mag <= sinkPar3) arg = 1.0;
          }
          else if(sinkShapeType == 2){	//exponential
            eps = sinkPar3;
            eps = pow(eps, sinkPar4);
            epsfactor = sinkPar5;
            if(epsfactor <= 0.0) epsfactor = 1.0;
            double magPow = pow(mag, sinkPar4);
            arg = exp(-magPow/(eps*epsfactor));
          }
          else if(sinkShapeType == 1){	//polynomial, compact support
            eps = sinkPar3;
            double pwrM = sinkPar4;
            double pwrN = sinkPar5;
            arg = 1.0 - pow((mag/eps),pwrM);
            arg = pow(arg, pwrN);
            if (mag >= eps){
              arg = 0.0;
            }
          }
          //the part that depends on sink rate profile
          double rate = 0.0;
          rate = sinkPar1*thePlanets[pi].omega; //constant
          if(sinkRateType == 2){  // sink rate depends on viscous time at sink radius
            double xnu[3];
            double rnu = thePlanets[pi].r + sinkPar3;
            xnu[0] = rnu;
            xnu[1] = thePlanets[pi].phi;
            xnu[2] = 0.0;
            double nu = get_nu(xnu, prim);
            double tvisc = 2*(mag2 + thePlanets[pi].eps*thePlanets[pi].eps)/(nu*3.0);
            rate = rate/tvisc;
          }
          else if(sinkRateType == 1){  // sink rate depends on viscous time locally
            double nu = get_nu(x, prim);
            double tvisc = 2*(mag2 + thePlanets[pi].eps*thePlanets[pi].eps)/(nu*3.0);
            rate = rate/tvisc;
          }

          surfdiff = rho*rate*arg;

          double delta = fmin(sinkPar2, 1.0);
          delta = fmax(0.0, sinkPar2);
          double rp, omp, vxp, vyp, vxg, vyg, vxr, vyr, vp_p, vp_r, vxn, vyn, cphi, sphi, vg_r, vg_p;
          rp = thePlanets[pi].r;
          omp = thePlanets[pi].omega;
          vp_p = rp*omp;
          vp_r = thePlanets[pi].vr;
          vxp = vp_r*cosp - vp_p*sinp;
          vyp = vp_r*sinp + vp_p*cosp;

          double vxg1 = vr*cosg - vp*sing;
          double vyg1 = vr*sing + vp*cosg;

          vxr = vxg1 - vxp;
          vyr = vyg1 - vyp;
          cphi = dx/mag;
          sphi = dy/mag;

          double acc_factor = dV*dt*surfdiff;

          double vpr = cphi*vyr - sphi*vxr;
          thePlanets[pi].Ls += (1.0-delta)*mag*vpr*acc_factor;
          vxn = (cphi*cphi + (1.0-delta)*sphi*sphi)*vxr + delta*sphi*cphi*vyr;
          vyn = delta*cphi*sphi*vxr + (sphi*sphi + (1.0-delta)*cphi*cphi)*vyr;

          vxg = vxn + vxp;
          vyg = vyn + vyp;
          vg_r =  vxg*cosg + vyg*sing;
          vg_p = -vxg*sing + vyg*cosg;

          thePlanets[pi].accL += acc_factor*(px*vyg - py*vxg + dx*vyp - dy*vxp);

          cons[DDD] -= acc_factor;
          cons[SRR] -= vg_r*acc_factor;
          cons[LLL] -= r*vg_p*acc_factor;
          cons[SZZ] -= vz*acc_factor;
          double v2 = vxg*vxg + vyg*vyg;
          cons[TAU] -= acc_factor*(specenth + 0.5*v2 - 0.5*(pow(vxg-vxg1, 2.0) + pow(vyg-vyg1, 2.0)) );

          thePlanets[pi].dM += acc_factor;

          double dxp = thePlanets[1-pi].r*cos(thePlanets[1-pi].phi) - px;
          double dyp = thePlanets[1-pi].r*sin(thePlanets[1-pi].phi) - py;

          //grouping these together avoids some cancellations and preserves precision when time-averaging
          double accEpart = -(dx*dxp + dy*dyp)*thePlanets[1-pi].M*pow(rbin, -3.0);
          accEpart -= 0.5*(vxp*vxp + vyp*vyp);
          accEpart -= thePlanets[1-pi].M/rbin;
          accEpart += vxg*vxp + vyg*vyp;
          thePlanets[pi].accE += accEpart*acc_factor;



          //not actually a sink, just torque accounting
          double fr,fp,fz;
          planetaryForce( thePlanets + pi, r, phi,  0.0, &fr, &fp, &fz, 0);
          thePlanets[pi].gravL -= rho*dV*dt*(px*(fr*sing + fp*cosg) - py*(fr*cosg - fp*sing));
          thePlanets[pi].gravE -= rho*dV*dt*(vxp*(fr*cosg - fp*sing) + vyp*(fr*sing + fp*cosg));
      }
    }
}




void cooling(double *prim, double *cons, double *xp, double *xm, double dV, double dt )
{
  if(coolType == COOL_BETA || coolType == COOL_BETA_RELAX)
  {
    //Beta-cooling
    double press = prim[PPP];
    double rho = prim[RHO];
    double gm1 = gamma_law-1.0;
    double beta = coolPar1;
    double x[3];
    get_centroid_arr(xp, xm, x);
    double enTarget = get_cs2(x)/(gamma_law*gm1);
    double enCurrent = press/(rho*gm1);
    double omtot = 0.0;

    if(coolType == COOL_BETA_RELAX || enCurrent > enTarget)
    {
      double r = x[0];
      double phi = x[1];
      double z = x[2];

      double cosg = cos(phi);
      double sing = sin(phi);
      double gx = r*cosg;
      double gy = r*sing;

      int pi;
      double fr,fp,fz, cosp, sinp, px, py, dx, dy, mag;
      for (pi=0; pi<Npl; pi++)
      {
        cosp = cos(thePlanets[pi].phi);
        sinp = sin(thePlanets[pi].phi);
        px = thePlanets[pi].r*cosp;
        py = thePlanets[pi].r*sinp;

        dx = gx-px;
        dy = gy-py;
        mag = dx*dx + dy*dy;
        mag = sqrt(mag);
        planetaryForce( thePlanets + pi, r, phi, z, &fr, &fp, &fz, 1);
        omtot += sqrt(fr*fr + fp*fp + fz*fz)/mag;
      }
      omtot = sqrt(omtot);

      ////direct implementation of source term
      //cons[TAU] += rho*(enCurrent - enTarget)*dt*dV*beta*omtot;

      //integrate source term over timestep
      double Tm1 = expm1(-dt*omtot/beta);
      cons[TAU] += rho*dV*( enCurrent - enTarget)*Tm1;	//N.B. Tm1 is in [-1 and 0]
    }
  }
}


void damping(double *prim, double *cons, double *xp, double *xm, double dV, double dt )
{
  if (DAMP_INNER + DAMP_OUTER + DAMP_LOWER + DAMP_UPPER > 0){
    double x[3];
    get_centroid_arr(xp, xm, x);

    double omtot = 1.0;
    if (DAMP_INNER==2 || DAMP_OUTER==2 || DAMP_LOWER==2 || DAMP_UPPER==2){
      double r = x[0];
      double phi = x[1];
      double z = x[2];

      double cosg = cos(phi);
      double sing = sin(phi);
      double gx = r*cosg;
      double gy = r*sing;

      int pi;
      omtot = 0.0;
      double fr,fp,fz, cosp, sinp, px, py, dx, dy, mag;
      for (pi=0; pi<Npl; pi++)
      {
        cosp = cos(thePlanets[pi].phi);
        sinp = sin(thePlanets[pi].phi);
        px = thePlanets[pi].r*cosp;
        py = thePlanets[pi].r*sinp;

        dx = gx-px;
        dy = gy-py;
        mag = dx*dx + dy*dy;
        mag = sqrt(mag);
        planetaryForce( thePlanets + pi, r, phi, z, &fr, &fp, &fz, 1);
        omtot += sqrt(fr*fr + fp*fp + fz*fz)/mag;
      }
      omtot = sqrt(omtot);
    }
    double ratetot = 0.0;
    double count = 0.0;
    double dampTime;
    double dampFactor;
    double dampLen;
    double theta;
    if (DAMP_INNER > 0){
      dampTime = dampTimeInner;
      if (DAMP_INNER == 2) dampTime = dampTime/omtot;
      dampLen = dampLenInner;
      theta = (x[0]-rmin)/dampLen;
      if (theta > 1.0) dampFactor = 0.0;
      else dampFactor = 1.0 - pow(1.0 - pow(theta,2.0), 2.0);
      ratetot = (count*ratetot + dampFactor/dampTime)/(1.0+count);
      count = count + 1.0;
    }
    if (DAMP_OUTER > 0){
      dampTime = dampTimeOuter;
      if (DAMP_OUTER == 2) dampTime = dampTime/omtot;
      dampLen = dampLenOuter;
      theta = (rmax-x[0])/dampLen;
      if (theta > 1.0) dampFactor = 0.0;
      else dampFactor = pow(1.0 - pow(theta,2.0), 2.0);
      ratetot = (count*ratetot + dampFactor/dampTime)/(1.0+count);
      count = count + 1.0;
    }
    if (twoD == 0){
      if (DAMP_UPPER > 0){
        dampTime = dampTimeUpper;
        if (DAMP_UPPER == 2) dampTime = dampTime/omtot;
        dampLen = dampLenUpper;
        theta = (zmax-x[2])/dampLen;
        if (theta > 1.0) dampFactor = 0.0;
        else dampFactor = pow(1.0 - pow(theta,2.0), 2.0);
        ratetot = (count*ratetot + dampFactor/dampTime)/(1.0+count);
        count = count + 1.0;
      }
      if (DAMP_LOWER > 0){
        dampTime = dampTimeLower;
        if (DAMP_LOWER == 2) dampTime = dampTime/omtot;
        dampLen = dampLenLower;
        theta = (zmax-x[2])/dampLen;
        if (theta > 1.0) dampFactor = 0.0;
        else dampFactor = 1.0-pow(1.0 - pow(theta,2.0), 2.0);
        ratetot = (count*ratetot + dampFactor/dampTime)/(1.0+count);
        count = count + 1.0;
      }
    }
    dampFactor = expm1(-dt*dampFactor);
    //dampFactor = -1.0*dt*dampFactor;
    double prims0[NUM_Q];
    double cons0[NUM_Q];
    double cons1[NUM_Q];
    initial(prims0, x);
    prim2cons(prims0, cons0, x, dV);
    prim2cons(prim, cons1, x, dV);
    cons[DDD] += (cons1[DDD] - cons0[DDD])*dampFactor;
    cons[SRR] += (cons1[SRR] - cons0[SRR])*dampFactor;
    cons[LLL] += (cons1[LLL] - cons0[LLL])*dampFactor;
    cons[SZZ] += (cons1[SZZ] - cons0[SZZ])*dampFactor;
    cons[TAU] += (cons1[TAU] - cons0[TAU])*dampFactor;
  }
}
