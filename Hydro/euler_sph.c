
#include "../paul.h"
#include "../hydro.h"
#include "../geometry.h"
#include "../omega.h"

static double gamma_law = 0.0; 
static double RHO_FLOOR = 0.0; 
static double PRE_FLOOR = 0.0; 
static int include_viscosity = 0;
static int isothermal = 0;
static int polar_sources_r = 0;
static int polar_sources_th = 0;

void setHydroParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
   isothermal = theDomain->theParList.isothermal_flag;
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
   include_viscosity = theDomain->theParList.visc_flag;
   if(theDomain->theParList.NoBC_Rmin == 1)
       polar_sources_r = 1;
   if(theDomain->theParList.NoBC_Zmin == 1
        || theDomain->theParList.NoBC_Zmax == 1)
       polar_sources_th = 1;
}

int set_B_flag(void){
   return(0);
}

double get_omega( const double * prim , const double * x ){
   return( prim[UPP] );
}


void prim2cons( const double * prim , double * cons , const double * x ,
                double dV, const double *xp, const double *xm ){

   double r = x[0];
   double sinth = sin(x[2]);
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r*sinth;
   double vt  = prim[UZZ]*r;
   double om  = get_om( x );
   double vp_off = vp - om*r*sinth;

   double v2  = vr*vr + vp_off*vp_off + vt*vt;

   double rhoe = Pp/(gamma_law - 1.);

   cons[DDD] = rho*dV;
   cons[TAU] = (.5*rho*v2 + rhoe )*dV;
   cons[SRR] = rho*vr*dV;
   cons[LLL] = r*sinth*rho*vp*dV;
   cons[SZZ] = r*rho*vt*dV;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      cons[q] = prim[q]*cons[DDD];
   }
}

void getUstar( const double * prim , double * Ustar , const double * x , double Sk , double Ss , const double * n , const double * Bpack ){

   double r = x[0];
   double sinth = sin(x[2]);
   double rho = prim[RHO];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r*sinth;
   double vt  = prim[UZZ]*r;
   double Pp  = prim[PPP];

   double om = get_om( x );
   double vp_off = vp - om*r*sinth;
   double v2 = vr*vr+vp_off*vp_off+vt*vt;

   double vn = vr*n[0] + vp*n[1] + vt*n[2];

   double vn_off = vn - om*r*sinth*n[1];
   double Ss_off = Ss + vn_off - vn;

   double rhoe = Pp/(gamma_law - 1.);

   double rhostar = rho*(Sk - vn)/(Sk - Ss);
   double Pstar = Pp*(Ss - vn)/(Sk - Ss);
   double Us = rhoe*(Sk - vn)/(Sk - Ss);

   Ustar[DDD] = rhostar;
   Ustar[SRR] =         rhostar*( vr + (Ss-vn)*n[0] );
   Ustar[LLL] = r*sinth*rhostar*( vp + (Ss-vn)*n[1] );
   Ustar[SZZ] =       r*rhostar*( vt + (Ss-vn)*n[2] );
   Ustar[TAU] = .5*rhostar*v2 + Us + rhostar*Ss_off*(Ss - vn) + Pstar;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      Ustar[q] = prim[q]*Ustar[DDD];
   }

}

void cons2prim( const double * cons , double * prim , const double * x , double dV, const double *xp, const double *xm ){
   
   double r = x[0];
   double sinth = sin(x[2]);
   double rho = cons[DDD]/dV;
   if( rho < RHO_FLOOR )   rho = RHO_FLOOR;
   double Sr  = cons[SRR]/dV;
   double Sp  = cons[LLL]/(dV*r*sinth);
   double St  = cons[SZZ]/(dV*r);
   double E   = cons[TAU]/dV;
   double om  = get_om( x );
   
   double vr = Sr/rho;
   double vp = Sp/rho;
   double vp_off = vp - om*r*sinth;
   double vt = St/rho;

   double KE = .5*( Sr*vr + rho*vp_off*vp_off + St*vt );
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
   prim[UPP] = vp/(r*sinth);
   prim[UZZ] = vt/r;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      prim[q] = cons[q]/cons[DDD];
   }

   //printf("c2p @ (%.6lg %.6lg %.6lg)\n", x[0], x[1], x[2]);
   //printf("    cons: %.6lg %.6lg %.6lg %.6lg %.6lg\n",
   //             cons[0], cons[1], cons[2], cons[3], cons[4]);
   //printf("    prim: %.6lg %.6lg %.6lg %.6lg %.6lg\n",
   //             prim[0], prim[1], prim[2], prim[3], prim[4]);

}

void flux( const double * prim , double * flux , const double * x , const double * n, const double *xp, const double *xm ){
  
   double r = x[0];
   double sinth = sin(x[2]);
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r*sinth;
   double vt  = prim[UZZ]*r;
   double om  = get_om( x );

   double vn = vr*n[0] + vp*n[1] + vt*n[2];
   double wn = om*r*sinth*n[1];
   double vp_off = vp - om*r*sinth;

   double rhoe = Pp/(gamma_law - 1.);
   double v2 = vr*vr + vp_off*vp_off + vt*vt;

   flux[DDD] = rho*vn;
   flux[SRR] = rho*vr*vn + Pp*n[0];
   flux[LLL] = r*sinth*(rho*vp*vn + Pp*n[1]);
   flux[SZZ] = r*(rho*vt*vn + Pp*n[2]);
   flux[TAU] = ( .5*rho*v2 + rhoe + Pp )*vn - Pp*wn;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      flux[q] = prim[q]*flux[DDD];
   }
   /*
   for(q=0; q<NUM_Q; q++)
       if(flux[q] != flux[q])
           printf("Whoa flux (%.0lg %.0lg %.0lg) NaN @ (%.6lg %.6lg %.6lg)\n",
                    n[0], n[1], n[2], x[0], x[1], x[2]);
                    */
   
}

void source( const double * prim , double * cons , const double * xp , const double * xm , double dVdt )
{
   
   double x[3];
   get_centroid_arr(xp, xm, x);
   double r = x[0];
   double th = x[2];
   double sinth = sin(th);
   double costh = cos(th);

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double ur = prim[URR];
   double up = prim[UPP];
   double ut = prim[UZZ];

   //double x[3] = {r, xm[1]+0.5*dphi, th};

   //Polar_Sources are the result of integrating the centripetal source term
   //in a cartesian frame, assuming rho and omega are constant. This leads to
   //better behaviour at r=0.
   //
   //TODO: IMPLEMENT THIS
   //
   //The naive source term (polar_sources==0), on the other hand, may exactly
   //cancel with source terms from other physics in the code. Gravitational
   //source terms obey polar_sources.
   //
   double centrifugal_r = rho*r*sinth*sinth*up*up + rho*r*ut*ut;
   double centrifugal_th = rho*r*r*sinth*costh*up*up;

   if(polar_sources_r || polar_sources_th)
   {
      double adjust[3];
      geom_polar_vec_adjust(xp, xm, adjust);
      if(polar_sources_r)
          centrifugal_r *= adjust[0];
      if(polar_sources_th)
          centrifugal_th *= adjust[2];
   }

   double r1_2 = 0.5*(xp[0]+xm[0]);
   double r2_3 = (xp[0]*xp[0] + xp[0]*xm[0] + xm[0]*xm[0]) / 3.0;
   double r3_4 = (3.0*xp[0]*r2_3 + xm[0]*xm[0]*xm[0]) / 4.0;
   double sinth1_2 = sin(0.5*(xp[2]+xm[2]));
   double costh1_2 = cos(0.5*(xp[2]+xm[2]));

   double press_bal_r, press_bal_th;
   press_bal_r = 2*Pp*r1_2/r2_3;
   press_bal_th = Pp*costh1_2/sinth1_2 * (r3_4*r1_2/(r2_3*r2_3));

   cons[SRR] += dVdt*( centrifugal_r + press_bal_r );
   cons[SZZ] += dVdt*( centrifugal_th + press_bal_th );


   // Om for energy frame
   double om  = get_om( x );
   double om_r = get_om1( x );
   double om_t = get_om2( x );

   cons[TAU] += dVdt*rho*( r*sinth*(sinth*ur+r*costh*ut) * om*om
                            - r*r*sinth*sinth*(up-om) * (ur*om_r + ut*om_t));
}

void visc_flux(const double * prim, const double * gradr, const double * gradp,
               const double * gradt, double * flux,
               const double * x, const double * n)
{
   double r = x[0];
   double th = x[2];
   double sinth = sin(th);
   double costh = cos(th);
   double nu = get_nu(x, prim);

   double rho = prim[RHO];
   double vr  = prim[URR];
   double om  = prim[UPP];
   double om_off = om - get_om(x);
   double vt  = prim[UZZ];

   //Divergence of v divided by number of spatial dimensions (3)
   double divV_o_d = (gradr[URR] + gradp[UPP] + gradt[UZZ] + 2*vr/r
                      + costh*vt/sinth) / 3.0;

   // Covariant components of shear tensor.
   double srr = gradr[URR] - divV_o_d;
   double spp = r*r*sinth*sinth*(gradp[UPP] + vr/r + costh*vt/sinth - divV_o_d);
   double stt = r*(r*gradt[UZZ] + vr - r*divV_o_d);
   double srp = 0.5*(r*r*sinth*sinth*gradr[UPP] + gradp[URR]);
   double srt = 0.5*(r*r*gradr[UZZ] + gradt[URR]);
   double spt = 0.5*r*r*(gradp[UZZ] + sinth*sinth*gradt[UPP]);

   // Covariant components of shear normal to face, shear_{ij} * n^{j}.
   // Given n is in orthonormal basis, 1/r factor corrects to coordinate basis
   double nc[3] = {n[0], n[1]/(r*sinth), n[2]/r};
   double srn = srr*nc[0] + srp*nc[1] + srt*nc[2];
   double spn = srp*nc[0] + spp*nc[1] + spt*nc[2];
   double stn = srt*nc[0] + spt*nc[1] + stt*nc[2];

   flux[SRR] = -2 * nu * rho * srn;
   flux[LLL] = -2 * nu * rho * spn;
   flux[SZZ] = -2 * nu * rho * stn;
   flux[TAU] = -2 * nu * rho * ( vr*srn + om_off*spn + vt*stn);
}

void visc_source(const double * prim, const double * gradr, const double *gradp,
                 const double * gradt, double * cons, const double *xp,
                 const double *xm, double dVdt)
{
   double x[3];
   get_centroid_arr(xp, xm, x);
   double r = x[0];
   double th = x[2];
   double sinth = sin(th);
   double costh = cos(th);
   double nu = get_nu(x, prim);

   double rho = prim[RHO];
   double vr  = prim[URR];
   double vt  = prim[UZZ];

   //Divergence of v divided by number of spatial dimensions (3)
   double divV_o_d = (gradr[URR] + gradp[UPP] + gradt[UZZ] + 2*vr/r
                      + costh*vt/sinth) / 3.0;

   // Relevant contravariant components of shear tensor.
   double spp = (r*gradp[UPP] + vr + r*costh*vt/sinth - r*divV_o_d)
                    / (r*r*r*sinth*sinth);
   double stt = (r*gradt[UZZ] + vr - r*divV_o_d) / (r*r*r);

   cons[SRR] += (-2 * rho * nu * (r * stt + r*sinth*sinth * spp)) * dVdt;
   cons[SZZ] += (-2 * rho * nu * (r*r*sinth*costh * spp)) * dVdt;

   // Mixed ^r_phi and ^theta_phi components of shear tensor
   double srp = 0.5*(r*r*sinth*sinth*gradr[UPP] + gradp[URR]);
   double spt = 0.5*(gradp[UZZ] + sinth*sinth*gradt[UPP]);
   double om_r = get_om1( x );
   double om_t = get_om2( x );

   cons[TAU] += (2 * rho * nu * (srp * om_r + spt * om_t)) * dVdt;
}

void flux_to_E( const double * Flux , const double * Ustr , const double * x , double * E1_riemann , double * B1_riemann , double * E2_riemann , double * B2_riemann , int dim ){

   //Silence is Golden.

}

void vel( const double * prim1 , const double * prim2 , double * Sl , double * Sr , double * Ss , const double * n , const double * x , double * Bpack ){

   double r = x[0];
   double sinth = sin(x[2]);
   double P1   = prim1[PPP];
   double rho1 = prim1[RHO];
   double vn1  = prim1[URR]*n[0] + prim1[UPP]*n[1]*r*sinth + prim1[UZZ]*n[2]*r;

   double cs1 = sqrt(gamma_law*P1/rho1);

   double P2   = prim2[PPP];
   double rho2 = prim2[RHO];
   double vn2  = prim2[URR]*n[0] + prim2[UPP]*n[1]*r*sinth + prim2[UZZ]*n[2]*r;

   double cs2 = sqrt(gamma_law*P2/rho2);

   *Ss = ( P2 - P1 + rho1*vn1*(-cs1) - rho2*vn2*cs2 )/( rho1*(-cs1) - rho2*cs2 );

   *Sr =  cs1 + vn1;
   *Sl = -cs1 + vn1;

   if( *Sr <  cs2 + vn2 ) *Sr =  cs2 + vn2;
   if( *Sl > -cs2 + vn2 ) *Sl = -cs2 + vn2;

}


double mindt(const double * prim , double w , const double * xp , const double * xm ){

   double x[3];
   get_centroid_arr(xp, xm, x);
   double r = x[0];
   double sinth = sin(x[2]);
   double Pp  = prim[PPP];
   double rho = prim[RHO];
   double vp  = (prim[UPP]-w)*r*sinth;
   double vr  = prim[URR];
   double vt  = prim[UZZ]*r;
   double cs  = sqrt(gamma_law*Pp/rho);

   double maxvr = cs + fabs(vr);
   double maxvp = cs + fabs(vp);
   double maxvt = cs + fabs(vt);

   double dtr = get_dL(xp,xm,1)/maxvr;
   double dtp = get_dL(xp,xm,0)/maxvp;
   double dtth= get_dL(xp,xm,2)/maxvt;
   
   double dt = dtr;
   if( dt > dtp ) dt = dtp;
   if( dt > dtth ) dt = dtth;

   if(include_viscosity)
   {
       double dL0 = get_dL(xp,xm,0);
       double dL1 = get_dL(xp,xm,1);
       double dL2 = get_dL(xp,xm,2);
       
       double dx = dL0;
       if( dx>dL1 ) dx = dL1;
       if( dx>dL2 ) dx = dL2;

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
