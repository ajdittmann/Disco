
#include "../paul.h"
#include "../geometry.h"
#include "../planet.h"

#define SNAP_NUM_Q 19
#define SNAP_NUM_R 30
#define SNAP_MAX_R 1.0

static double gamma_law = 0.0;

void setDiagParams( struct domain * theDomain )
{
   gamma_law = theDomain->theParList.Adiabatic_Index;
}

int num_diagnostics(void)
{
    return(15);
}

int num_snapshot_rz(void)
{
    return SNAP_NUM_Q;
}

int num_snapshot_arr(void)
{
    return 2*(SNAP_NUM_Q+1)*SNAP_NUM_R;
}

/* Diagnostics for Circumbinary Disks */

void get_diagnostics(const double *x, const double *prim, double *Qrz, 
                        struct domain * theDomain )
{
    double rpz[3];
    get_rpz(x, rpz);
    double r = rpz[0];
    double phi = rpz[1];
    double z = rpz[2];

    double rho = prim[RHO];
    double Pp = prim[PPP];

    double V[3] = {prim[URR], prim[UPP], prim[UZZ]};
    get_vec_covariant(x, V, V); // V is now in the orthonormal basis
    double Vrpz[3];
    get_vec_rpz(x, V, Vrpz);
    double vr = Vrpz[0];
    double vp = Vrpz[1];
    double vz = Vrpz[2];


    Qrz[0] = rho;           // Sigma
    Qrz[1] = rho*V[0];        // Mdot
    Qrz[2] = rho*fabs(V[0]);  // Mdot |vr|
    Qrz[3] = rho * r*vp;    // angular momentum
    Qrz[4] = rho*r*vp*V[0];	// (advective) angular momentum flux
    Qrz[5] = rho*r*vp*fabs(V[0]); 
    Qrz[6] = Pp;

    double cosp = cos(phi);
    double sinp = sin(phi);
    double xyz[3] = {r*cosp, r*sinp, z};

    double vx = vr*cosp - vp*sinp;
    double vy = vr*sinp + vp*cosp;
    double v2 =  vr*vr + vp*vp + vz*vz;
    double rdv = r * vr + z * vz;
    double sinTh = r / sqrt(r*r + z*z);
    Qrz[7] = rho * ((r*v2 - sinTh)*cosp - rdv*vx);	//e_x
    Qrz[8] = rho * ((r*v2 - sinTh)*sinp - rdv*vy);	//e_y

    double Fxyz[3];
    double Fp;

    planetaryForce( theDomain->thePlanets + 0, xyz, Fxyz);
    Fp = cosp * Fxyz[1] - sinp * Fxyz[0];
    Qrz[9] = rho * r * Fp;				//Torque density from pl 0
    planetaryForce( theDomain->thePlanets + 1, xyz, Fxyz);
    Fp = cosp * Fxyz[1] - sinp * Fxyz[0];
    Qrz[10] = rho * r * Fp;				//Torque density from pl 1

    double cos2p = (cosp - sinp) * (cosp + sinp);
    double sin2p = 2*sinp*cosp;
    Qrz[11]  = rho * r * cosp;
    Qrz[12]  = rho * r * sinp;
    Qrz[13]  = rho * r*r * cos2p;
    Qrz[14] = rho * r*r * sin2p;
}


void get_snapshot_rz(const double *x, const double *prim, double *Qrz, 
                        struct domain * theDomain )
{
    double rpz[3];
    get_rpz(x, rpz);
    double r = rpz[0];
    double phi = rpz[1];
    double z = rpz[2];

    double cosp = cos(phi);
    double sinp = sin(phi);
    double xyz[3] = {r*cosp, r*sinp, z};

    double rho = prim[RHO];
    double V[3] = {prim[URR], prim[UPP], prim[UZZ]};
    double Pp = prim[PPP];

    get_vec_covariant(x, V, V); // V is now in the orthonormal basis
    double Vrpz[3];
    get_vec_rpz(x, V, Vrpz);
    double vr = Vrpz[0];
    double vp = Vrpz[1];
    double vz = Vrpz[2];

    double vx = vr*cosp - vp*sinp;
    double vy = vr*sinp + vp*cosp;
    double v2 =  vr*vr + vp*vp + vz*vz;
    double rdv = r * vr + z * vz;
    double sinTh = r / sqrt(r*r + z*z);

    double F1xyz[3];
    double F2xyz[3];

    planetaryForce( theDomain->thePlanets + 0, xyz, F1xyz);
    double Fp1 = cosp * F1xyz[1] - sinp * F1xyz[0];
    planetaryForce( theDomain->thePlanets + 1, xyz, F2xyz);
    double Fp2 = cosp * F2xyz[1] - sinp * F2xyz[0];

    double cos2p = (cosp - sinp) * (cosp + sinp);
    double sin2p = 2*sinp*cosp;
    
    Qrz[0] = rho;           // Sigma
    Qrz[1] = rho*V[0];        // Mdot
    Qrz[2] = rho * r*vp;    // angular momentum
    Qrz[3] = rho*r*vp*V[0];	// (advective) angular momentum flux
    Qrz[4] = Pp;
    Qrz[5] = rho * ((r*v2 - sinTh)*cosp - rdv*vx);	//e_x
    Qrz[6] = rho * ((r*v2 - sinTh)*sinp - rdv*vy);	//e_y
    Qrz[7] = rho * r * Fp1;				//Torque density from pl 0
    Qrz[8] = rho * F1xyz[0];			//Force density from pl 0
    Qrz[9] = rho * F1xyz[1];			//Force density from pl 0
    Qrz[10] = rho * F1xyz[2];			//Force density from pl 0
    Qrz[11] = rho * r * Fp2;				//Torque density from pl 1
    Qrz[12] = rho * F2xyz[0];			//Force density from pl 1
    Qrz[13] = rho * F2xyz[1];			//Force density from pl 1
    Qrz[14] = rho * F2xyz[2];			//Force density from pl 1
    Qrz[15]  = rho * cosp;
    Qrz[16] = rho * sinp;
    Qrz[17] = rho * cos2p;
    Qrz[18] = rho * sin2p;
}

void get_snapshot_arr(const double *x, const double *prim, double *Qarr, 
                        struct domain * theDomain )
{
    memset(Qarr, 0, 2*(SNAP_NUM_Q+1)*SNAP_NUM_R * sizeof(double));

    int p;
    for(p=0; p < 2; p++)
    {
        double rpz[3];
        get_rpz(x, rpz);
        double r = rpz[0];
        double phi = rpz[1];
        double z = rpz[2];

        double cosp = cos(phi);
        double sinp = sin(phi);
        double xyz[3] = {r*cosp, r*sinp, z};

        struct planet *pl = theDomain->thePlanets + p;

        double dx = xyz[0] - pl->xyz[0];
        double dy = xyz[1] - pl->xyz[1];

        double dr = sqrt(dx*dx + dy*dy);

        if(dr > SNAP_MAX_R)
            continue;

        int j = (int) ((dr * SNAP_NUM_R)/SNAP_MAX_R);
        if(j < 0 || j >= SNAP_NUM_R)
            continue;

        int idx = (p*SNAP_NUM_R + j) * (SNAP_NUM_Q+1);

        //Now we know we are near-enough to the minidisk to compute things.

        double rho = prim[RHO];
        double V[3] = {prim[URR], prim[UPP], prim[UZZ]};
        double Pp = prim[PPP];

        get_vec_covariant(x, V, V); // V is now in the orthonormal basis
        double Vrpz[3];
        get_vec_rpz(x, V, Vrpz);
        double vx = Vrpz[0]*cosp - Vrpz[1]*sinp;
        double vy = Vrpz[0]*sinp + Vrpz[1]*cosp;
        double vz = Vrpz[2];

        double cosp_pl = pl->xyz[0] / pl->r;
        double sinp_pl = pl->xyz[1] / pl->r;

        double vxp = cosp_pl * pl->vr - sinp_pl * pl->r * pl->omega;
        double vyp = sinp_pl * pl->vr + cosp_pl * pl->r * pl->omega;

        // Re-base all coordinates to be relative to the planet
        
        r = dr;
        cosp = dx / dr;
        sinp = dy / dr;

        vx -= vxp;
        vy -= vyp;

        // Now can calculate things in the planet frame

        double vr = vx*cosp + vy*sinp;
        double vp = -vx*sinp + vy*cosp;

        double v2 =  vx*vx + vy*vy + vz*vz;
        double rdv = r * vr + z * vz;
        double sinTh = r / sqrt(r*r + z*z);

        double F1xyz[3];
        double F2xyz[3];

        planetaryForce( theDomain->thePlanets + 0, xyz, F1xyz);
        double Fp1 = cosp * F1xyz[1] - sinp * F1xyz[0];
        planetaryForce( theDomain->thePlanets + 1, xyz, F2xyz);
        double Fp2 = cosp * F2xyz[1] - sinp * F2xyz[0];

        double cos2p = (cosp - sinp) * (cosp + sinp);
        double sin2p = 2*sinp*cosp;
        double mpl = pl->M;

        Qarr[idx+0] = 1.0;
        Qarr[idx+1] = rho;           // Sigma
        Qarr[idx+2] = rho*vr;        // Mdot
        Qarr[idx+3] = rho * r*vp;    // angular momentum
        Qarr[idx+4] = rho*r*vp*vr;	// (advective) angular momentum flux
        Qarr[idx+5] = Pp;
        Qarr[idx+6] = rho * ((r*v2 - mpl*sinTh)*cosp - rdv*vx) / mpl;	//e_x
        Qarr[idx+7] = rho * ((r*v2 - mpl*sinTh)*sinp - rdv*vy) / mpl;	//e_y
        Qarr[idx+8] = rho * r * Fp1;				//Torque density from pl 0
        Qarr[idx+9] = rho * F1xyz[0];			//Force density from pl 0
        Qarr[idx+10] = rho * F1xyz[1];			//Force density from pl 0
        Qarr[idx+11] = rho * F1xyz[2];			//Force density from pl 0
        Qarr[idx+12] = rho * r * Fp2;				//Torque density from pl 1
        Qarr[idx+13] = rho * F2xyz[0];			//Force density from pl 1
        Qarr[idx+14] = rho * F2xyz[1];			//Force density from pl 1
        Qarr[idx+15] = rho * F2xyz[2];			//Force density from pl 1
        Qarr[idx+16]  = rho * cosp;
        Qarr[idx+17] = rho * sinp;
        Qarr[idx+18] = rho * cos2p;
        Qarr[idx+19] = rho * sin2p;
    }
}
