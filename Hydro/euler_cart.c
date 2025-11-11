
#include "../paul.h"
#include "../hydro.h"
#include "../geometry.h"
#include "../omega.h"

static double gamma_law = 0.0; 
static double RHO_FLOOR = 0.0; 
static double PRE_FLOOR = 0.0; 
static int include_viscosity = 0;
static int isothermal = 0;

void setHydroParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
   isothermal = theDomain->theParList.isothermal_flag;
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
   include_viscosity = theDomain->theParList.visc_flag;
}

int set_B_flag(void){
   return(0);
}

double get_omega( const double * prim , const double * x ){
   return( prim[UPP] );
}


void prim2cons( const double * prim , double * cons , const double * x , double dV, const double *xp, const double *xm ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[URR];
   double vp  = prim[UPP];
   double vz  = prim[UZZ];
   double om  = get_om( x );
   double vp_off = vp - om;

   double v2  = vr*vr + vp_off*vp_off + vz*vz;

   double rhoe = Pp/(gamma_law - 1.);

   cons[DDD] = rho*dV;
   cons[TAU] = (.5*rho*v2 + rhoe )*dV;
   cons[SRR] = rho*vr*dV;
   cons[LLL] = rho*vp*dV;
   cons[SZZ] = rho*vz*dV;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      cons[q] = prim[q]*cons[DDD];
   }
}

void getUstar( const double * prim , double * Ustar , const double * x , double Sk , double Ss , const double * n , const double * Bpack ){

   double rho = prim[RHO];
   double vr  = prim[URR];
   double vp  = prim[UPP];
   double vz  = prim[UZZ];
   double Pp  = prim[PPP];

   double om = get_om( x );
   double vp_off = vp - om;
   double v2 = vr*vr+vp_off*vp_off+vz*vz;

   double vn = vr*n[0] + vp*n[1] + vz*n[2];

   double vn_off = vn - om*n[1];
   double Ss_off = Ss + vn_off - vn;

   double rhoe = Pp/(gamma_law - 1.);

   double rhostar = rho*(Sk - vn)/(Sk - Ss);
   double Pstar = Pp*(Ss - vn)/(Sk - Ss);
   double Us = rhoe*(Sk - vn)/(Sk - Ss);

   Ustar[DDD] = rhostar;
   Ustar[SRR] = rhostar*( vr + (Ss-vn)*n[0] );
   Ustar[LLL] = rhostar*( vp + (Ss-vn)*n[1] );
   Ustar[SZZ] = rhostar*( vz + (Ss-vn)*n[2] );
   Ustar[TAU] = .5*rhostar*v2 + Us + rhostar*Ss_off*(Ss - vn) + Pstar;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      Ustar[q] = prim[q]*Ustar[DDD];
   }

}

void cons2prim( const double * cons , double * prim , const double * x , double dV, const double *xp, const double *xm ){
   
   double rho = cons[DDD]/dV;
   if( rho < RHO_FLOOR )   rho = RHO_FLOOR;
   double Sr  = cons[SRR]/dV;
   double Sp  = cons[LLL]/dV;
   double Sz  = cons[SZZ]/dV;
   double E   = cons[TAU]/dV;
   double om  = get_om( x );
   
   double vr = Sr/rho;
   double vp = Sp/rho;
   double vp_off = vp - om;
   double vz = Sz/rho;

   double KE = .5*( Sr*vr + rho*vp_off*vp_off + Sz*vz );
   double rhoe = E-KE;
   double Pp = (gamma_law - 1.)*rhoe;

   if( Pp  < PRE_FLOOR*rho ) Pp = PRE_FLOOR*rho;
   if( isothermal ){
      double cs2 = get_cs2( x );
      Pp = cs2*rho/gamma_law;
   }

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr;
   prim[UPP] = vp;
   prim[UZZ] = vz;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      prim[q] = cons[q]/cons[DDD];
   }

}

void flux( const double * prim , double * flux , const double * x , const double * n, const double *xp, const double *xm ){
   
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[URR];
   double vp  = prim[UPP];
   double vz  = prim[UZZ];
   double om  = get_om( x );

   double vn = vr*n[0] + vp*n[1] + vz*n[2];
   double wn = om*n[1];
   double vp_off = vp - om;

   double rhoe = Pp/(gamma_law - 1.);
   double v2 = vr*vr + vp_off*vp_off + vz*vz;

   flux[DDD] = rho*vn;
   flux[SRR] = rho*vr*vn + Pp*n[0];
   flux[LLL] = rho*vp*vn + Pp*n[1];
   flux[SZZ] = rho*vz*vn + Pp*n[2];
   flux[TAU] = ( .5*rho*v2 + rhoe + Pp )*vn - Pp*wn;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      flux[q] = prim[q]*flux[DDD];
   }
   
}

void source( const double * prim , double * cons , const double * xp , const double * xm , double dVdt ){}

void visc_flux(const double * prim, const double * gradx, const double * grady,
               const double * gradz, double * flux,
               const double * x, const double * n)
{
   double nu = get_nu(x, prim);

   double rho = prim[RHO];
   double vx  = prim[URR];
   double vy  = prim[UPP];
   double vy_off = vy - get_om(x);
   double vz  = prim[UZZ];

   //Divergence of v divided by number of spatial dimensions (3)
   double divV_o_d = (gradx[URR] + grady[UPP] + gradz[UZZ]) / 3.0;

   // Covariant components of shear tensor.
   double sxx = gradx[URR] - divV_o_d;
   double syy = grady[UPP] - divV_o_d;
   double szz = gradz[UZZ] - divV_o_d;
   double sxy = 0.5*(gradx[UPP] + grady[URR]);
   double sxz = 0.5*(gradx[UZZ] + gradz[URR]);
   double syz = 0.5*(grady[UZZ] + gradz[UPP]);

   // Covariant components of shear normal to face, shear_{ij} * n^{j}.
   // Given n is in orthonormal basis, 1/r factor corrects to coordinate basis
   double sxn = sxx*n[0] + sxy*n[1] + sxz*n[2];
   double syn = sxy*n[0] + syy*n[1] + syz*n[2];
   double szn = sxz*n[0] + syz*n[1] + szz*n[2];

   flux[SRR] = -2 * nu * rho * sxn;
   flux[LLL] = -2 * nu * rho * syn;
   flux[SZZ] = -2 * nu * rho * szn;
   flux[TAU] = -2 * nu * rho * ( vx*sxn + vy_off*syn + vz*szn);
}

void visc_source(const double * prim, const double * gradr, const double *gradp,
                 const double * gradt, double * cons, const double *xp,
                 const double *xm, double dVdt){}

void flux_to_E( const double * Flux , const double * Ustr , const double * x , double * E1_riemann , double * B1_riemann , double * E2_riemann , double * B2_riemann , int dim ){

   //Silence is Golden.

}

void vel( const double * prim1 , const double * prim2 , double * Sl , double * Sr , double * Ss , const double * n , const double * x , double * Bpack ){

   double P1   = prim1[PPP];
   double rho1 = prim1[RHO];
   double vn1  = prim1[URR]*n[0] + prim1[UPP]*n[1] + prim1[UZZ]*n[2];

   double cs1 = sqrt(gamma_law*P1/rho1);

   double P2   = prim2[PPP];
   double rho2 = prim2[RHO];
   double vn2  = prim2[URR]*n[0] + prim2[UPP]*n[1] + prim2[UZZ]*n[2];

   double cs2 = sqrt(gamma_law*P2/rho2);

   *Ss = ( P2 - P1 + rho1*vn1*(-cs1) - rho2*vn2*cs2 )/( rho1*(-cs1) - rho2*cs2 );

   *Sr =  cs1 + vn1;
   *Sl = -cs1 + vn1;

   if( *Sr <  cs2 + vn2 ) *Sr =  cs2 + vn2;
   if( *Sl > -cs2 + vn2 ) *Sl = -cs2 + vn2;

}


double mindt(const double * prim , double w , const double * xp , const double * xm ){

   double Pp  = prim[PPP];
   double rho = prim[RHO];
   double vp  = (prim[UPP]-w);
   double vr  = prim[URR];
   double vz  = prim[UZZ];
   double cs  = sqrt(gamma_law*Pp/rho);

   double maxvr = cs + fabs(vr);
   double maxvp = cs + fabs(vp);
   double maxvz = cs + fabs(vz);

   double dtr = get_dL(xp,xm,1)/maxvr;
   double dtp = get_dL(xp,xm,0)/maxvp;
   double dtz = get_dL(xp,xm,2)/maxvz;
   
   double dt = dtr;
   if( dt > dtp ) dt = dtp;
   if( dt > dtz ) dt = dtz;

   if(include_viscosity)
   {
       double dL0 = get_dL(xp,xm,0);
       double dL1 = get_dL(xp,xm,1);
       double dL2 = get_dL(xp,xm,2);
       
       double dx = dL0;
       if( dx>dL1 ) dx = dL1;
       if( dx>dL2 ) dx = dL2;

       double x[3];
       get_centroid_arr(xp, xm, x);
       double nu = get_nu(x, prim);

       double dt_visc = 0.5*dx*dx/nu;
       if( dt > dt_visc )
           dt = dt_visc;
   }

   return( dt );
}

void reflect_prims(double * prim, const double * x, int dim)
{
    //dim == 0: r, dim == 1: p, dim == 2: z
    if(dim == 0)
        prim[URR] = -prim[URR];
    else if(dim == 1)
        prim[UPP] = -prim[UPP];
    else if(dim == 2)
        prim[UZZ] = -prim[UZZ];
}

double bfield_scale_factor(double x, int dim)
{
    // Returns the factor used to scale B_cons.
    // x is coordinate location in direction dim.
    // dim == 0: r, dim == 1: p, dim == 2: z
    
    return 1.0;
}

double getCartInterpWeight(const double *x)
{
    return 0.0;
}
