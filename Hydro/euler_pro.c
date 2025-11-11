
#include "../paul.h"
#include "../hydro.h"
#include "../geometry.h"
#include "../omega.h"

static double gamma_law = 0.0; 
static double RHO_FLOOR = 0.0; 
static double PRE_FLOOR = 0.0; 
static int include_viscosity = 0;
static int isothermal = 0;
static int polar_sources_mu = 0;
static int polar_sources_nu = 0;

void setHydroParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
   isothermal = theDomain->theParList.isothermal_flag;
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
   include_viscosity = theDomain->theParList.visc_flag;
   if(theDomain->theParList.NoBC_Rmin == 1)
       polar_sources_mu = 1;
   if(theDomain->theParList.NoBC_Zmin == 1
        || theDomain->theParList.NoBC_Zmax == 1)
       polar_sources_nu = 1;
}

int set_B_flag(void){
   return(0);
}

double get_omega( const double * prim , const double * x ){
   return( prim[UPP] );
}


void prim2cons( const double * prim , double * cons , const double * x ,
                double dV, const double *xp, const double *xm )
{
    double sh_mu = sinh(x[0]);
    double s_nu = sin(x[2]);
    double ht = sqrt(sh_mu*sh_mu + s_nu*s_nu);
    double hp = sh_mu * s_nu;
   
    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double vm  = prim[URR] * ht;
    double vp  = prim[UPP] * hp;
    double vn  = prim[UZZ] * ht;
    double om  = get_om( x );
    double vp_off = vp - om*hp;

    double v2  = vm*vm + vp_off*vp_off + vn*vn;

    double rhoe = Pp/(gamma_law - 1.);

    cons[DDD] = rho*dV;
    cons[TAU] = (.5*rho*v2 + rhoe )*dV;
    cons[SRR] = ht*rho*vm*dV;
    cons[LLL] = hp*rho*vp*dV;
    cons[SZZ] = ht*rho*vn*dV;

    int q;
    for( q=NUM_C ; q<NUM_Q ; ++q ){
        cons[q] = prim[q]*cons[DDD];
    }
}

void getUstar( const double * prim , double * Ustar , const double * x , double Sk , double Ss , const double * n , const double * Bpack )
{
    double sh_mu = sinh(x[0]);
    double s_nu = sin(x[2]);
    double ht = sqrt(sh_mu*sh_mu + s_nu*s_nu);
    double hp = sh_mu * s_nu;
   
    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double vm  = prim[URR] * ht;
    double vp  = prim[UPP] * hp;
    double vn  = prim[UZZ] * ht;

    double om = get_om( x );
    double vp_off = vp - om*hp;
    double v2 = vm*vm + vp_off*vp_off + vn*vn;

    double vN = vm*n[0] + vp*n[1] + vn*n[2];

    double vN_off = vN - om*hp*n[1];
    double Ss_off = Ss + vN_off - vN;

    double rhoe = Pp/(gamma_law - 1.);

    double rhostar = rho*(Sk - vN)/(Sk - Ss);
    double Pstar = Pp*(Ss - vN)/(Sk - Ss);
    double Us = rhoe*(Sk - vN)/(Sk - Ss);

    Ustar[DDD] = rhostar;
    Ustar[SRR] = ht * rhostar*( vm + (Ss-vN)*n[0] );
    Ustar[LLL] = hp * rhostar*( vp + (Ss-vN)*n[1] );
    Ustar[SZZ] = ht * rhostar*( vn + (Ss-vN)*n[2] );
    Ustar[TAU] = .5*rhostar*v2 + Us + rhostar*Ss_off*(Ss - vN) + Pstar;

    int q;
    for( q=NUM_C ; q<NUM_Q ; ++q ){
        Ustar[q] = prim[q]*Ustar[DDD];
    }
}

void cons2prim(const double *cons, double *prim, const double *x, double dV,
               const double *xp, const double *xm)
{   
    double sh_mu = sinh(x[0]);
    double s_nu = sin(x[2]);
    double ht = sqrt(sh_mu*sh_mu + s_nu*s_nu);
    double hp = sh_mu * s_nu;
    
    double rho = cons[DDD]/dV;
    if( rho < RHO_FLOOR )
        rho = RHO_FLOOR;
   
    double Sm  = cons[SRR]/(dV*ht);
    double Sp  = cons[LLL]/(dV*hp);
    double Sn  = cons[SZZ]/(dV*ht);
    double E   = cons[TAU]/dV;
    double om  = get_om( x );
   
    double vm = Sm/rho;
    double vp = Sp/rho;
    double vn = Sn/rho;
    double vp_off = vp - om * hp;

    double KE = .5*( Sm*vm + rho*vp_off*vp_off + Sn*vn );
    double rhoe = E-KE;
    double Pp = (gamma_law - 1.)*rhoe;

    if( Pp  < PRE_FLOOR*rho )
        Pp = PRE_FLOOR*rho;

    if( isothermal ){
        double cs2 = get_cs2( x );
        Pp = cs2*rho/gamma_law;
    }

    prim[RHO] = rho;
    prim[PPP] = Pp;
    prim[URR] = vm / ht;
    prim[UPP] = vp / hp;
    prim[UZZ] = vn / ht;

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

void flux(const double *prim, double *flux, const double *x, const double *n,
          const double *xp, const double *xm)
{  
    double sh_mu = sinh(x[0]);
    double s_nu = sin(x[2]);
    double ht = sqrt(sh_mu*sh_mu + s_nu*s_nu);
    double hp = sh_mu * s_nu;
   
    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double vm  = prim[URR] * ht;
    double vp  = prim[UPP] * hp;
    double vn  = prim[UZZ] * ht;
    double om  = get_om( x );
    double vp_off = vp - om*hp;

    double vN = vm*n[0] + vp*n[1] + vn*n[2];
    double wN = om*hp*n[1];

    double rhoe = Pp/(gamma_law - 1.);
    double v2 = vm*vm + vp_off*vp_off + vn*vn;

    flux[DDD] = rho*vN;
    flux[SRR] = ht * (rho*vm*vN + Pp*n[0]);
    flux[LLL] = hp * (rho*vp*vN + Pp*n[1]);
    flux[SZZ] = ht * (rho*vn*vN + Pp*n[2]);
    flux[TAU] = ( .5*rho*v2 + rhoe + Pp )*vN - Pp*wN;

    int q;
    for(q=NUM_C; q<NUM_Q; ++q)
        flux[q] = prim[q]*flux[DDD];
}

void source(const double *prim, double *cons, 
            const double *xp, const double *xm, double dVdt)
{
    double x[3];
    get_centroid_arr(xp, xm, x);

    double sh_mu = sinh(x[0]);
    double ch_mu = cosh(x[0]);
    double s_nu = sin(x[2]);
    double c_nu = cos(x[2]);
    //double ht = sqrt(sh_mu*sh_mu + s_nu*s_nu);
    double hp = sh_mu * s_nu;
   
    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double um = prim[URR];
    double up = prim[UPP];
    double un = prim[UZZ];

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

    double centrifugal_mu = rho * sh_mu*ch_mu
            * (um*um + un*un + s_nu*s_nu * up*up);
    double centrifugal_nu = rho * s_nu*c_nu
            * (um*um + un*un + sh_mu*sh_mu * up*up);

    if(polar_sources_mu || polar_sources_nu)
    {
        double adjust[3];
        geom_polar_vec_adjust(xp, xm, adjust);
        if(polar_sources_mu)
            centrifugal_mu *= adjust[0];
        if(polar_sources_nu)
            centrifugal_nu *= adjust[2];
    }

    double sh_mup = sinh(xp[0]);
    double sh_mum = sinh(xm[0]);
    double s_nup = sin(xp[0]);
    double s_num = sin(xm[0]);
    double dA_mup = sqrt(sh_mup*sh_mup + s_nu*s_nu) * sh_mup;
    double dA_mum = sqrt(sh_mum*sh_mum + s_nu*s_nu) * sh_mum;
    double dV_mu = (sh_mu*sh_mu + s_nu*s_nu) * sh_mu * (xp[0]-xm[0]);
    double dA_nup = sqrt(sh_mu*sh_mu + s_nup*s_nup) * s_nup;
    double dA_num = sqrt(sh_mu*sh_mu + s_num*s_num) * s_num;
    double dV_nu = (sh_mu*sh_mu + s_nu*s_nu) * s_nu * (xp[2]-xm[2]);

    double press_bal_mu = Pp * (dA_mup - dA_mum) / dV_mu;
    double press_bal_nu = Pp * (dA_nup - dA_num) / dV_nu;

    cons[SRR] += dVdt*( centrifugal_mu + press_bal_mu );
    cons[SZZ] += dVdt*( centrifugal_nu + press_bal_nu );

    // Om for energy frame
    double om  = get_om( x );
    double om_m = get_om1( x );
    double om_n = get_om2( x );

    cons[TAU] += dVdt*rho*(
            sh_mu*s_nu * (ch_mu*s_nu * um + sh_mu*c_nu * un) * om*om
                - hp*hp*(up-om) * (um*om_m + un*om_n));
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

void vel(const double *prim1, const double *prim2, double  *Sl, double *Sr,
         double *Ss, const double *n, const double *x, double *Bpack)
{
    double sh_mu = sinh(x[0]);
    double s_nu = sin(x[2]);
    double ht = sqrt(sh_mu*sh_mu + s_nu*s_nu);
    double hp = sh_mu * s_nu;
   
    double rho1 = prim1[RHO];
    double P1  = prim1[PPP];
    double vn1 = ht*prim1[URR]*n[0] + hp*prim1[UPP]*n[1] + ht*prim1[UZZ]*n[2];

    double cs1 = sqrt(gamma_law*P1/rho1);

    double P2   = prim2[PPP];
    double rho2 = prim2[RHO];
    double vn2 = ht*prim2[URR]*n[0] + hp*prim2[UPP]*n[1] + ht*prim2[UZZ]*n[2];

    double cs2 = sqrt(gamma_law*P2/rho2);

    *Ss = ( P2 - P1 + rho1*vn1*(-cs1) - rho2*vn2*cs2 )
            / ( rho1*(-cs1) - rho2*cs2 );

    *Sr =  cs1 + vn1;
    *Sl = -cs1 + vn1;

    if( *Sr <  cs2 + vn2 )
        *Sr =  cs2 + vn2;
    if( *Sl > -cs2 + vn2 )
        *Sl = -cs2 + vn2;
}


double mindt(const double *prim, double w, const double *xp, const double *xm)
{
    double x[3];
    get_centroid_arr(xp, xm, x);

    double sh_mu = sinh(x[0]);
    double s_nu = sin(x[2]);
    double ht = sqrt(sh_mu*sh_mu + s_nu*s_nu);
    double hp = sh_mu * s_nu;

    double Pp  = prim[PPP];
    double rho = prim[RHO];
    double vm  = prim[URR] * ht;
    double vp  = (prim[UPP]-w) * hp;
    double vn  = prim[UZZ] * ht;
    double cs  = sqrt(gamma_law*Pp/rho);

    double maxvm = cs + fabs(vm);
    double maxvp = cs + fabs(vp);
    double maxvn = cs + fabs(vn);

    double dtm = get_dL(xp,xm,1) / maxvm;
    double dtp = get_dL(xp,xm,0) / maxvp;
    double dtn = get_dL(xp,xm,2) / maxvn;
   
    double dt = dtp;
    if( dt > dtm ) dt = dtm;
    if( dt > dtn ) dt = dtn;

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
