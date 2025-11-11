#include "../paul.h"
#include "../hydro.h"
#include "../geometry.h"
#include "../omega.h"

static double gamma_law = 0.0; 
static double mach = 1.0; 
static double RHO_FLOOR = 0.0; 
static double PRE_FLOOR = 0.0; 
static int include_viscosity = 0;
static int isothermal = 0;
static int polar_sources = 0;
#if ENABLE_CART_INTERP
static int Cartesian_Interp = 0;
static double Cartesian_Interp_R0 = 0;
#endif

void setHydroParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
   isothermal = theDomain->theParList.isothermal_flag;
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
   mach = theDomain->theParList.Disk_Mach;
   include_viscosity = theDomain->theParList.visc_flag;
#if ENABLE_CART_INTERP
   Cartesian_Interp = theDomain->theParList.Cartesian_Interp;
   Cartesian_Interp_R0 = theDomain->theParList.Cartesian_Interp_R0;
#endif
   if(theDomain->theParList.NoBC_Rmin == 1)
       polar_sources = 1;
}

double get_nu( const double *, const double *);

int set_B_flag(void){
   return(0);
}

double get_omega( const double * prim , const double * x ){
    UNUSED(x);
   return( prim[UPP] );
}

void prim2cons(const double * prim, double * cons, const double * x,
               double dV, const double * xp, const double *xm) {

   double r = x[0];
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r;
   double vz  = prim[UZZ];
   double om  = get_om( x );
   double vp_off = vp - om*r;

   double v2  = vr*vr + vp_off*vp_off + vz*vz;

   double rhoe = Pp/(gamma_law - 1.);

   double cart_adjust = 1.0;

#if ENABLE_CART_INTERP
   if(xp != NULL && xm != NULL && Cartesian_Interp)
   {
      double dphi = xp[1] - xm[1];
      cart_adjust = 2*sin(0.5*dphi) / dphi;
   }
#endif

   cons[DDD] = rho*dV;
   cons[TAU] = (.5*rho*v2 + rhoe )*dV;
   cons[SRR] = rho*vr*dV * cart_adjust;
   cons[LLL] = r*rho*vp*dV * cart_adjust;
   cons[SZZ] = rho*vz*dV;

   if(0 && NUM_N > 0)
   {
       // EXPERIMENTAL. Ideal-gas specific entropy (NOT VALID FOR RAD-PRESSURE GAS)
       double s = log(Pp * pow(rho, -gamma_law)) / (gamma_law - 1.0);
       cons[NUM_C] = rho * s * dV;
   }

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      cons[q] = prim[q]*cons[DDD];
   }
}

void getUstar( const double * prim , double * Ustar , const double * x ,
              double Sk , double Ss , const double * n , const double * Bpack ){

    UNUSED(Bpack);

   double r = x[0];
   double rho = prim[RHO];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r;
   double vz  = prim[UZZ];
   double Pp  = prim[PPP];

   double om = get_om( x );
   double vp_off = vp - om*r;
   double v2 = vr*vr+vp_off*vp_off+vz*vz;

   double vn = vr*n[0] + vp*n[1] + vz*n[2];

   double vn_off = vn - om*r*n[1];
   double Ss_off = Ss + vn_off - vn;

   double rhoe = Pp/(gamma_law - 1.);

   double rhostar = rho*(Sk - vn)/(Sk - Ss);
   double Pstar = Pp*(Ss - vn)/(Sk - Ss);
   double Us = rhoe*(Sk - vn)/(Sk - Ss);

   Ustar[DDD] = rhostar;
   Ustar[SRR] =   rhostar*( vr + (Ss-vn)*n[0] );
   Ustar[LLL] = r*rhostar*( vp + (Ss-vn)*n[1] );
   Ustar[SZZ] =   rhostar*( vz + (Ss-vn)*n[2] );
   Ustar[TAU] = .5*rhostar*v2 + Us + rhostar*Ss_off*(Ss - vn) + Pstar;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      Ustar[q] = prim[q]*Ustar[DDD];
   }

}

void cons2prim( const double * cons , double * prim , const double * x ,
                double dV, const double *xp, const double *xm){
   
   double r = x[0];
   double cart_adjust = 1.0;
#if ENABLE_CART_INTERP
   if(xp != NULL && xm != NULL && Cartesian_Interp)
   {
      double dphi = xp[1] - xm[1];
      cart_adjust = 2*sin(0.5*dphi) / dphi;
   }
#endif

   double rho = cons[DDD]/dV;
   if( rho < RHO_FLOOR )   rho = RHO_FLOOR;
   double Sr  = cons[SRR]/(dV * cart_adjust);
   double Sp  = cons[LLL]/(dV * r * cart_adjust);
   double Sz  = cons[SZZ]/dV;
   double E   = cons[TAU]/dV;
   double om  = get_om( x );
   
   double vr = Sr/rho;
   double vp = Sp/rho;
   double vp_off = vp - om*r;
   double vz = Sz/rho;

   double KE = .5*( Sr*vr + rho*vp_off*vp_off + Sz*vz );
   double rhoe = E-KE;
   double Pp = (gamma_law - 1.)*rhoe;

   if(0 && (Pp <= 0.0 && NUM_N > 0))
   {
       //EXPERIMENTAL ENTROPY FIX
       double s = cons[NUM_C] / cons[DDD];
       Pp = exp((gamma_law-1.0)*s) * pow(rho, gamma_law);
   }
   else if(Pp < PRE_FLOOR*rho)
   {
       Pp = PRE_FLOOR*rho;
   }

   if( isothermal ){
      double cs2 = get_cs2( x );
      Pp = cs2*rho/gamma_law;
   }

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr;
   prim[UPP] = vp/r;
   prim[UZZ] = vz;

   if(0 && NUM_N > 0)
   {
       //EXPERIMENTAL entropy fix
       double s = log(Pp * pow(rho, -gamma_law)) / (gamma_law - 1.0);
       prim[NUM_C] = s;
   }

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      prim[q] = cons[q]/cons[DDD];
   }

}

void flux(const double * prim, double * flux, const double * x,
          const double * n, const double *xp, const double *xm)
{
    double r = x[0];
    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double vr  = prim[URR];
    double vp  = prim[UPP]*r;
    double vz  = prim[UZZ];
    double om  = get_om( x );

    double vn = vr*n[0] + vp*n[1] + vz*n[2];
    double wn = om*r*n[1];
    double vp_off = vp - om*r;

    double rhoe = Pp/(gamma_law - 1.);
    double v2 = vr*vr + vp_off*vp_off + vz*vz;

    int q;

#if ENABLE_CART_INTERP
    if(xp != NULL && xm != NULL && Cartesian_Interp)
    {
        double r1 = 0.5*(xp[0] + xm[0]);
        double dphi = xp[1] - xm[1];
        double scalar_fac = dphi > 0.0 ? 2 * sin(0.5*dphi) / dphi : 1.0;
        double vector_fac = dphi > 0.0 ? sin(dphi) / dphi : 1.0;
        double rad_fac_1 = 0.5*(1 + vector_fac);
        double rad_fac_0 = 0.5*(1 - vector_fac);
        double scal_fac_n = scalar_fac * n[0] + n[1] + n[2];
        
        flux[DDD] = rho*vn * scal_fac_n;
        flux[SRR] = (rad_fac_1 * rho*vr*vr + rad_fac_0 * rho*vp*vp + Pp) * n[0]
                    + (rho * vr * vp)*n[1]
                    + (rho*vr*vz * scalar_fac) * n[2];
        flux[LLL] = r*rho*vr*vp * vector_fac * n[0]
                    + r1*(rho*vp*vp + Pp) * n[1]
                    + r*rho*vp*vz * vector_fac * n[2];
        flux[SZZ] = rho*vz*vr * vector_fac * n[0]
                    + rho*vz*vp * n[1]
                    + (rho*vz*vz + Pp) * n[2];
        flux[TAU] = ((0.5*rho*v2 + rhoe + Pp)*vn - Pp*wn) * scal_fac_n;
    }
    else
    {
#endif
        flux[DDD] = rho*vn;
        flux[SRR] = rho*vr*vn + Pp*n[0];
        flux[LLL] = r*(rho*vp*vn + Pp*n[1]);
        flux[SZZ] = rho*vz*vn + Pp*n[2];
        flux[TAU] = ( .5*rho*v2 + rhoe + Pp )*vn - Pp*wn;

#if ENABLE_CART_INTERP
    }
#endif
    for( q=NUM_C ; q<NUM_Q ; ++q )
        flux[q] = prim[q]*flux[DDD];
}

void source( const double * prim , double * cons , const double * xp , const double * xm , double dVdt)
{
   double rp = xp[0];
   double rm = xm[0];
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double r = get_centroid(rp, rm, 1);
   double z = get_centroid(xp[2], xm[2], 2);
   double r_1  = .5*(rp+rm);
   double vr  = prim[URR];
   double omega = prim[UPP];
   double vz  = prim[UZZ];

   double x[3] = {r, 0.5*(xm[1]+xp[1]), z};

   //Polar_Sources are the result of integrating the centripetal source term
   //in a cartesian frame, assuming rho and omega are constant. This leads to
   //better behaviour at r=0.
   //
   //The naive source term (polar_sources==0), on the other hand, may exactly
   //cancel with source terms from other physics in the code. Gravitational
   //source terms obey polar_sources.
   //
   double centrifugal = rho*omega*omega*r;
   
   if(polar_sources)
   {
      double adjust[3];
      geom_polar_vec_adjust(xp, xm, adjust);
      centrifugal *= adjust[0];
   }
#if ENABLE_CART_INTERP
   if(Cartesian_Interp)
   {
        double dphi = xp[1] - xm[1];
        //double scalar_fac = 2 * sin(0.5*dphi) / dphi;
        double vector_fac = sin(dphi) / dphi;
        double rad_fac_1 = 0.5*(1 + vector_fac);
        double rad_fac_0 = 0.5*(1 - vector_fac);
        double vp = r * omega;
        centrifugal = (rho*vr*vr * rad_fac_0 + rho*vp*vp * rad_fac_1) / r_1;
   }
#endif

   double press_bal   = Pp/r_1;

   // Geometric source term
   cons[SRR] += dVdt*( centrifugal + press_bal );

   // Om for energy frame
   double om  = get_om( x );
   double om_r = get_om1( x );
   double om_z = get_om2( x );
   
   // energy-adjustment source
   cons[TAU] += dVdt*rho*( r*vr*om*om - r*r*(omega-om)*(vr*om_r + vz*om_z));
}

void visc_flux(const double * prim, const double * gradr, const double * gradp,
               const double * gradz, double * flux,
               const double * x, const double * n)
{
   double r = x[0];
   double nu = get_nu(x, prim);

   double rho = prim[RHO];
   double vr  = prim[URR];
   double om  = prim[UPP];
   double om_off = om - get_om(x);
   double vz  = prim[UZZ];

   //Divergence of v divided by number of spatial dimensions (3)
   double divV_o_d = (gradr[URR] + gradp[UPP] + gradz[UZZ] + vr/r) / 3.0;

   // Covariant components of shear tensor.
   double srr = gradr[URR] - divV_o_d;
   double spp = r*(r*gradp[UPP] + vr - r*divV_o_d);
   double szz = gradz[UZZ] - divV_o_d;
   double srp = 0.5*(r*r*gradr[UPP] + gradp[URR]);
   double srz = 0.5*(gradr[UZZ] + gradz[URR]);
   double spz = 0.5*(gradp[UZZ] + r*r*gradz[UPP]);

   // Covariant components of shear normal to face, shear_{ij} * n^{j}.
   // Given n is in orthonormal basis, 1/r factor corrects to coordinate basis
   double nc[3] = {n[0], n[1]/r, n[2]};
   double srn = srr*nc[0] + srp*nc[1] + srz*nc[2];
   double spn = srp*nc[0] + spp*nc[1] + spz*nc[2];
   double szn = srz*nc[0] + spz*nc[1] + szz*nc[2];

   flux[SRR] = -2 * nu * rho * srn;
   flux[LLL] = -2 * nu * rho * spn;
   flux[SZZ] = -2 * nu * rho * szn;
   flux[TAU] = -2 * nu * rho * ( vr*srn + om_off*spn + vz*szn);
}

void visc_source(const double * prim, const double * gradr, const double *gradp,
                 const double * gradz, double * cons, const double *xp,
                 const double *xm, double dVdt)
{
   double x[3];
   get_centroid_arr(xp, xm, x);
   double r = x[0];
   double nu = get_nu(x, prim);

   double rho = prim[RHO];
   double vr  = prim[URR];

   //Divergence of v divided by number of spatial dimensions (3)
   double divV_o_d = (gradr[URR] + gradp[UPP] + gradz[UZZ] + vr/r) / 3.0;

   // Contravariant -phi-phi component of shear tensor.
   double spp = (r*gradp[UPP] + vr - r*divV_o_d) / (r*r*r);
   
   cons[SRR] += (-2 * r * rho * nu * spp) * dVdt;

   // Mixed ^r_phi and ^z_phi components of shear tensor
   double srp = 0.5*(r*r*gradr[UPP] + gradp[URR]);
   double spz = 0.5*(gradp[UZZ] + r*r*gradz[UPP]);
   double om_r = get_om1( x );
   double om_z = get_om2( x );

   cons[TAU] += (2 * rho * nu * (srp * om_r + spz * om_z)) * dVdt;

   if(0 && NUM_N > 0)
   {
       //EXPERIMENTAL Viscous entropy source
       //
       // covariant _{ij} shear tensor components that we don't
       // already have
       double sRR = gradr[URR] - divV_o_d;
       double sZZ = gradz[UZZ] - divV_o_d;
       double sRZ = 0.5*(gradr[UZZ] + gradz[URR]);

       double s2 = sRR*sRR + spp*spp * (r*r*r*r) + sZZ*sZZ
                    + 2*(srp*srp + spz*spz)/(r*r) + 2*sRZ*sRZ;

       // IDEAL GAS TEMPERATURE. NOT APPROPRIATE FOR RAD-PRESSURE GAS
       double T = prim[PPP] / rho;

       cons[NUM_C] += 2*rho*nu * s2 / T * dVdt;
   }
}

void flux_to_E( const double * Flux , const double * Ustr , const double * x, 
                double * E1_riemann , double * B1_riemann , 
                double * E2_riemann , double * B2_riemann , int dim ){

   //Silence is Golden.

    UNUSED(Flux);
    UNUSED(Ustr);
    UNUSED(x);
    UNUSED(E1_riemann);
    UNUSED(B1_riemann);
    UNUSED(E2_riemann);
    UNUSED(B2_riemann);
    UNUSED(dim);

}

void vel( const double * prim1 , const double * prim2 , 
         double * Sl , double * Sr , double * Ss , 
         const double * n , const double * x , double * Bpack ){

    UNUSED(Bpack);

   double r = x[0];
   double P1   = prim1[PPP];
   double rho1 = prim1[RHO];
   double vn1  = prim1[URR]*n[0] + prim1[UPP]*n[1]*r + prim1[UZZ]*n[2];

   double cs1 = sqrt(gamma_law*P1/rho1);

   double P2   = prim2[PPP];
   double rho2 = prim2[RHO];
   double vn2  = prim2[URR]*n[0] + prim2[UPP]*n[1]*r + prim2[UZZ]*n[2];

   double cs2 = sqrt(gamma_law*P2/rho2);

   *Ss = ( P2 - P1 + rho1*vn1*(-cs1) - rho2*vn2*cs2 )/( rho1*(-cs1) - rho2*cs2 );

   *Sr =  cs1 + vn1;
   *Sl = -cs1 + vn1;

   if( *Sr <  cs2 + vn2 ) *Sr =  cs2 + vn2;
   if( *Sl > -cs2 + vn2 ) *Sl = -cs2 + vn2;

}

double mindt(const double * prim , double w ,
             const double * xp , const double * xm ){

   double x[3];
   get_centroid_arr(xp, xm, x);
   double r = x[0];
   double Pp  = prim[PPP];
   double rho = prim[RHO];
   double vp  = (prim[UPP]-w)*r;
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

       double nu = get_nu(x, prim);

       double dt_visc = 0.5*dx*dx/nu;
       if( dt > dt_visc )
           dt = dt_visc;
   }

   return( dt );

}
/*
double getReynolds( double * prim , double w , double * x , double dx ){
   double r = x[0];
   double nu = explicit_viscosity;
   double vr = prim[URR];
   double omega = prim[UPP];
   double vp = omega*r-w;
   double vz = prim[UZZ];
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double cs = sqrt(gamma_law*Pp/rho);
   
   double v = sqrt(vr*vr + vp*vp + vz*vz);
   double Re = (v+cs)*dx/nu;
   
   return(Re);
}
*/

void reflect_prims(double * prim, const double * x, int dim)
{
    UNUSED(x);

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
    
    UNUSED(x);
    UNUSED(dim);
    
    return 1.0;
}

#if ENABLE_CART_INTERP
double getCartInterpWeight(const double *x)
{
    if(x[0] >= Cartesian_Interp_R0)
        return 0.0;
    else if(x[0] >= 0.5 * Cartesian_Interp_R0)
        return (Cartesian_Interp_R0 - x[0]) / (0.5*Cartesian_Interp_R0);
    else
        return 1.0;
}
#endif
