
#include "../paul.h"
#include "../geometry.h"
#include "../planet.h"

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
    return 0;
}

int num_snapshot_arr(void)
{
    return 0;
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
    Qrz[2] = rho * r*vp;    // angular momentum
    Qrz[3] = rho*r*vp*V[0];	// (advective) angular momentum flux
    Qrz[4] = Pp;

    double cosp = cos(phi);
    double sinp = sin(phi);
    double xyz[3] = {r*cosp, r*sinp, z};

    double vx = vr*cosp - vp*sinp;
    double vy = vr*sinp + vp*cosp;
    double v2 =  vr*vr + vp*vp + vz*vz;
    double rdv = r * vr + z * vz;
    double sinTh = r / sqrt(r*r + z*z);
    Qrz[5] = rho * ((r*v2 - sinTh)*cosp - rdv*vx);	//e_x
    Qrz[6] = rho * ((r*v2 - sinTh)*sinp - rdv*vy);	//e_y

    double Fxyz[3];
    double Fp;

    planetaryForce( theDomain->thePlanets + 0, xyz, Fxyz);
    Fp = cosp * Fxyz[1] - sinp * Fxyz[0];
    Qrz[9] = rho * r * Fp;				//Torque density from pl 0
    Qrz[10] = rho * Fxyz[0];
    Qrz[11] = rho * Fxyz[1];
    Qrz[12] = rho * Fxyz[2];
    planetaryForce( theDomain->thePlanets + 1, xyz, Fxyz);
    Fp = cosp * Fxyz[1] - sinp * Fxyz[0];
    Qrz[13] = rho * r * Fp;				//Torque density from pl 1
    Qrz[14] = rho * Fxyz[0];
    Qrz[15] = rho * Fxyz[1];
    Qrz[16] = rho * Fxyz[2];

    double cos2p = (cosp - sinp) * (cosp + sinp);
    double sin2p = 2*sinp*cosp;
    Qrz[17]  = rho * cosp;
    Qrz[18]  = rho * sinp;
    Qrz[19]  = rho * cos2p;
    Qrz[20] = rho * sin2p;
}


void get_snapshot_rz(const double *x, const double *prim, double *Qrz, 
                        struct domain * theDomain )
{
    //Silence is golden
}

void get_snapshot_arr(const double *x, const double *prim, double *Qarr, 
                        struct domain * theDomain )
{
    //Silence is golden
}
