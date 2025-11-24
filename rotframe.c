#include "disco.h"
#include "omega.h"
#include "geometry.h"

static int om_flag = 0;
static double Omega0 = 0.0;
static double d = 0.0;

void setRotFrameParams( struct domain * theDomain ){
   om_flag = theDomain->theParList.RotFrame;
   Omega0  = theDomain->theParList.RotOmega;
   d       = theDomain->theParList.RotD;
}

void subtract_omega( double * prim ){
   if( om_flag ) prim[UPP] -= Omega0;
}

void omegaForce( double r , double phi , double vr , double omega , double * frpz ){

   //Omega0^2 vec R + 2 Omega0 x v

   double Rr = r-d*cos(phi);
   double Rp = d*sin(phi);

   frpz[0] = Omega0*Omega0*Rr + 2.*Omega0*r*omega;
   frpz[1] = Omega0*Omega0*Rp - 2.*Omega0*vr;
   frpz[2] = 0.0;

}

void omega_src( double * prim , double * cons , double * xp , double * xm , double dVdt ){

   if( om_flag ){
      /*
      double rho = prim[RHO];
      double vr  = prim[URR];
      double omega = prim[UPP];
   
      double r = get_centroid(xp[0], xm[0], 1);
      double vp  = r*omega;
      double dphi = get_dp(xp[1],xm[1]);
      double phi = xm[1] + 0.5*dphi;

      double Fr,Fp;
      omegaForce( r , phi , vr , omega , &Fr , &Fp );

      cons[SRR] += rho*Fr*dVdt;
      cons[LLL] += rho*Fp*r*dVdt;
      cons[TAU] += rho*( Fr*vr + Fp*vp )*dVdt;
      */

      double x[3], rpz[3];
      get_centroid_arr(xp, xm, x);
      get_rpz(x, rpz);

      double rho = prim[RHO];
      double v[3] = {prim[URR], prim[UPP], prim[UZZ]};
      double vrpz[3];
      get_vec_covariant(x, v, v); // v is now in orthonormal basis
      get_vec_rpz(x, v, vrpz);    // vrpz in orthonormal basis
      get_vec_contravariant(x, vrpz, vrpz);  //vrpz is contravariant

      double frpz[3], f[3];
      omegaForce(rpz[0], rpz[1], vrpz[0], vrpz[1], frpz);  //frpz is orthonormal
      get_vec_from_rpz(x, frpz, f);  //f is orthonormal

      //adjust v for energy-subtraction scheme.
      v[1] -= get_scale_factor(x, 0) * get_om(x);  //still orthonormal

      double vf = v[0]*f[0] + v[1]*f[1] + v[2]*f[2];  //both in same basis,
                                                      // easy dot product
      get_vec_covariant(x, f, f);  // f is now covariant

      cons[SRR] += rho * f[0] * dVdt;
      cons[LLL] += rho * f[1] * dVdt;
      cons[SZZ] += rho * f[2] * dVdt;
      cons[TAU] += rho * vf * dVdt;
   }

}

