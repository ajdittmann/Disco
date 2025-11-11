
#include "../paul.h"
#include "../hydro.h"
#include "../geometry.h"
#include "../omega.h"

static double gamma_law = 0.0; 
static double RHO_FLOOR = 0.0; 
static double PRE_FLOOR = 0.0; 
static int isothermal = 0;

void setHydroParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
   isothermal = theDomain->theParList.isothermal_flag;
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
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
   double v[3] = {prim[URR], prim[UPP], prim[UZZ]};

   double vxyz[3];
   get_vec_covariant(x, v, v);
   get_vec_xyz(x, v, vxyz);

   double v_off[3] = {prim[URR], prim[UPP]-get_om(x), prim[UZZ]};
   get_vec_covariant(x, v_off, v_off);
   double v2  = v_off[0]*v_off[0] + v_off[1]*v_off[1] + v_off[2]*v_off[2];

   double rhoe = Pp/(gamma_law - 1.);

   cons[DDD] = rho*dV;
   cons[TAU] = (.5*rho*v2 + rhoe )*dV;
   cons[SRR] = rho*vxyz[0]*dV;
   cons[LLL] = rho*vxyz[1]*dV;
   cons[SZZ] = rho*vxyz[2]*dV;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      cons[q] = prim[q]*cons[DDD];
   }
}

void getUstar( const double * prim , double * Ustar , const double * x , double Sk , double Ss , const double * n , const double * Bpack ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double v[3] = {prim[URR], prim[UPP], prim[UZZ]};
   double om = get_om(x);
   double v_off[3] = {prim[URR], prim[UPP]-om, prim[UZZ]};

   double vxyz[3], v_offxyz[3], nxyz[3];
   get_vec_covariant(x, v, v);
   get_vec_covariant(x, v_off, v_off);
   get_vec_xyz(x, v, vxyz);
   get_vec_xyz(x, v_off, v_offxyz);
   get_vec_xyz(x, n, nxyz);

   double v2  = v_offxyz[0]*v_offxyz[0] + v_offxyz[1]*v_offxyz[1]
                    + v_offxyz[2]*v_offxyz[2];

   double h[3] = {get_scale_factor(x,1), get_scale_factor(x,0),
                    get_scale_factor(x,2)};
   double vn = v[0]*n[0] + v[1]*n[1] + v[2]*n[2];

   double vn_off = vn - om*h[1]*n[1];
   double Ss_off = Ss + vn_off - vn;

   double rhoe = Pp/(gamma_law - 1.);

   double rhostar = rho*(Sk - vn)/(Sk - Ss);
   double Pstar = Pp*(Ss - vn)/(Sk - Ss);
   double Us = rhoe*(Sk - vn)/(Sk - Ss);

   Ustar[DDD] = rhostar;
   Ustar[SRR] = rhostar*( vxyz[0] + (Ss-vn)*nxyz[0] );
   Ustar[LLL] = rhostar*( vxyz[1] + (Ss-vn)*nxyz[1] );
   Ustar[SZZ] = rhostar*( vxyz[2] + (Ss-vn)*nxyz[2] );
   Ustar[TAU] = .5*rhostar*v2 + Us + rhostar*Ss_off*(Ss - vn) + Pstar;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      Ustar[q] = prim[q]*Ustar[DDD];
   }

}

void cons2prim( const double * cons , double * prim , const double * x , double dV, const double *xp, const double *xm ){
   
   double rho = cons[DDD]/dV;
   if( rho < RHO_FLOOR )   rho = RHO_FLOOR;
   double Sx  = cons[SRR]/dV;
   double Sy  = cons[LLL]/dV;
   double Sz  = cons[SZZ]/dV;
   double E   = cons[TAU]/dV;
   double om  = get_om( x );
   
   double vxyz[3] = {Sx/rho, Sy/rho, Sz/rho};

   double v[3];
   get_vec_from_xyz(x, vxyz, v);
   get_vec_contravariant(x, v, v);
   double v_off[3] = {v[0], v[1]-om, v[2]};
   get_vec_covariant(x, v_off, v_off);
   double v2 = v_off[0]*v_off[0] + v_off[1]*v_off[1] + v_off[2]*v_off[2];

   double KE = .5*rho*v2;
   double rhoe = E-KE;
   double Pp = (gamma_law - 1.)*rhoe;

   if( Pp  < PRE_FLOOR*rho ) Pp = PRE_FLOOR*rho;
   if( isothermal ){
      double cs2 = get_cs2( x );
      Pp = cs2*rho/gamma_law;
   }

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = v[0];
   prim[UPP] = v[1];
   prim[UZZ] = v[2];

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      prim[q] = cons[q]/cons[DDD];
   }

}


void flux( const double * prim , double * flux , const double * x , const double * n, const double *xp, const double *xm ){
   
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double om = get_om(x);

   double v[3] = {prim[URR], prim[UPP], prim[UZZ]};
   get_vec_covariant(x, v, v);
   double vxyz[3];
   get_vec_xyz(x, v, vxyz);
   
   double v_off[3] = {prim[URR], prim[UPP] - om, prim[UZZ]};
   get_vec_covariant(x, v_off, v_off);
   double v2 = v_off[0]*v_off[0] + v_off[1]*v_off[1] + v_off[2]*v_off[2];

   double nxyz[3];
   get_vec_xyz(x, n, nxyz);

   double vn = v[0]*n[0] + v[1]*n[1] + v[2]*n[2];
   double wn = get_scale_factor(x,0)*om*n[1];

   double rhoe = Pp/(gamma_law - 1.);

   flux[DDD] = rho*vn;
   flux[SRR] = rho*vxyz[0]*vn + Pp*nxyz[0];
   flux[LLL] = rho*vxyz[1]*vn + Pp*nxyz[1];
   flux[SZZ] = rho*vxyz[2]*vn + Pp*nxyz[2];
   flux[TAU] = ( .5*rho*v2 + rhoe + Pp )*vn - Pp*wn;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      flux[q] = prim[q]*flux[DDD];
   }
   
}


void source( const double * prim , double * cons , const double * xp , const double * xm , double dVdt ){}

void visc_flux(const double * prim, const double * gradr, const double * gradp,
               const double * gradz, double * flux,
               const double * x, const double * n){}
void visc_source(const double * prim, const double * gradr, const double *gradp,
                 const double * gradt, double * cons, const double *xp,
                 const double *xm, double dVdt){}

void flux_to_E( const double * Flux , const double * Ustr , const double * x , double * E1_riemann , double * B1_riemann , double * E2_riemann , double * B2_riemann , int dim ){

   //Silence is Golden.

}

void vel( const double * prim1 , const double * prim2 , double * Sl , double * Sr , double * Ss , const double * n , const double * x , double * Bpack ){

   double v1[3] = {prim1[URR], prim1[UPP], prim1[UZZ]};
   get_vec_covariant(x, v1, v1);

   double P1   = prim1[PPP];
   double rho1 = prim1[RHO];
   double vn1 = v1[0]*n[0] + v1[1]*n[1] + v1[2]*n[2];

   double cs1 = sqrt(gamma_law*P1/rho1);
   
   double v2[3] = {prim1[URR], prim1[UPP], prim1[UZZ]};
   get_vec_covariant(x, v2, v2);

   double P2   = prim2[PPP];
   double rho2 = prim2[RHO];
   double vn2 = v2[0]*n[0] + v2[1]*n[1] + v2[2]*n[2];

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
   double cs  = sqrt(gamma_law*Pp/rho);

   double x[3];
   get_centroid_arr(xp, xm, x);
   double v[3] = {prim[URR], prim[UPP], prim[UZZ]};
   v[1] -= w;

   double maxv1 = cs + fabs(v[0]);
   double maxv2 = cs + fabs(v[1]);
   double maxv3 = cs + fabs(v[2]);

   double dt1 = get_dL(xp,xm,1)/maxv1;
   double dt2 = get_dL(xp,xm,0)/maxv2;
   double dt3 = get_dL(xp,xm,2)/maxv3;
   
   double dt = dt1;
   if( dt > dt2 ) dt = dt2;
   if( dt > dt3 ) dt = dt3;

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
