
#include "../paul.h"
#include "../hydro.h"
#include "../geometry.h"
#include "../omega.h"

static double gamma_law = 0.0; 
static double RHO_FLOOR = 0.0; 
static double PRE_FLOOR = 0.0; 
static double explicit_viscosity = 0.0;
static int include_viscosity = 0;
static int isothermal = 0;
static int polar_sources = 0;

void setHydroParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
   isothermal = theDomain->theParList.isothermal_flag;
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
   explicit_viscosity = theDomain->theParList.viscosity;
   include_viscosity = theDomain->theParList.visc_flag;
   if(theDomain->NgRa == 0)
       polar_sources = 1;
}

int set_B_flag(void){
   return(1);
}

double get_omega( const double * prim , const double * x ){
   UNUSED(x);
   return( prim[UPP] );
}


void prim2cons( const double * prim , double * cons , const double * x , double dV, const double *xp, const double *xm){

   UNUSED(xp);
   UNUSED(xm);

   double r = x[0];
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r;
   double vz  = prim[UZZ];
   double om  = get_om( x );
   double vp_off = vp - om*r;

   double Br = prim[BRR];
   double Bp = prim[BPP];
   double Bz = prim[BZZ];

   double v2  = vr*vr + vp_off*vp_off + vz*vz;
   double B2  = Br*Br + Bp*Bp + Bz*Bz;

   double rhoe = Pp/(gamma_law - 1.);

   cons[DDD] = rho*dV;
   cons[TAU] = (.5*rho*v2 + rhoe + .5*B2 )*dV;
   cons[SRR] = rho*vr*dV;
   cons[LLL] = r*rho*vp*dV;
   cons[SZZ] = rho*vz*dV;

   cons[BRR] = Br*dV * bfield_scale_factor(r, 0);
   cons[BPP] = Bp*dV/r;
   cons[BZZ] = Bz*dV;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      cons[q] = prim[q]*cons[DDD];
   }
}

void getUstar( const double * prim , double * Ustar , const double * x , double Sk , double Ss , const double * n , const double * Bpack ){

   double r = x[0];

   double Bsn = Bpack[0];
   double Bsr = Bpack[1];
   double Bsp = Bpack[2];
   double Bsz = Bpack[3];
   double vBs = Bpack[4];

   double rho = prim[RHO];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r;
   double vz  = prim[UZZ];
   double Pp  = prim[PPP];

   double Br  = prim[BRR];
   double Bp  = prim[BPP];
   double Bz  = prim[BZZ];

   double v2 = vr*vr+vp*vp+vz*vz;
   double B2 = Br*Br+Bp*Bp+Bz*Bz;

   double vn = vr*n[0] + vp*n[1] + vz*n[2];
   double Bn = Br*n[0] + Bp*n[1] + Bz*n[2];
   double vB = vr*Br   + vp*Bp   + vz*Bz;

   double rhoe = Pp/(gamma_law-1.);

   double D  = rho; 
   double mr = rho*vr;
   double mp = rho*vp;
   double mz = rho*vz;
   double E  = .5*rho*v2 + rhoe + .5*B2;

   double Bs2 = Bsr*Bsr+Bsp*Bsp+Bsz*Bsz;
   double Ps  = rho*( Sk - vn )*( Ss - vn ) + (Pp+.5*B2-Bn*Bn) - .5*Bs2 + Bsn*Bsn;

   double Dstar = ( Sk - vn )*D/( Sk - Ss );
   double Msr   = ( ( Sk - vn )*mr + ( Br*Bn - Bsr*Bsn ) ) / ( Sk - Ss );
   double Msp   = ( ( Sk - vn )*mp + ( Bp*Bn - Bsp*Bsn ) ) / ( Sk - Ss );
   double Msz   = ( ( Sk - vn )*mz + ( Bz*Bn - Bsz*Bsn ) ) / ( Sk - Ss );
   double Estar = ( ( Sk - vn )*E + (Ps+.5*Bs2)*Ss - (Pp+.5*B2)*vn - vBs*Bsn + vB*Bn ) / ( Sk - Ss );

   double Msn = Dstar*Ss;
   double mn  = Msr*n[0] + Msp*n[1] + Msz*n[2];

   Msr += n[0]*( Msn - mn );
   Msp += n[1]*( Msn - mn );
   Msz += n[2]*( Msn - mn );

   Ustar[DDD] = Dstar;
   Ustar[SRR] = Msr;
   Ustar[LLL] = r*Msp;
   Ustar[SZZ] = Msz;
   Ustar[TAU] = Estar;

   Ustar[BRR] = Bsr * bfield_scale_factor(r, 0);
   Ustar[BPP] = Bsp/r;
   Ustar[BZZ] = Bsz;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      Ustar[q] = prim[q]*Ustar[DDD];
   }
}

void cons2prim( const double * cons , double * prim , const double * x , double dV, const double *xp, const double *xm ){
   UNUSED(xp);
   UNUSED(xm);

   double r = x[0];
   
   double rho = cons[DDD]/dV;
   if( rho < RHO_FLOOR )   rho = RHO_FLOOR;
   double Sr  = cons[SRR]/dV;
   double Sp  = cons[LLL]/dV/r;
   double Sz  = cons[SZZ]/dV;
   double E   = cons[TAU]/dV;
   double om  = get_om( x );
   
   double vr = Sr/rho;
   double vp = Sp/rho;
   double vp_off = vp - om*r;
   double vz = Sz/rho;

   double Br  = cons[BRR]/(dV * bfield_scale_factor(r, 0));
   double Bp  = cons[BPP]/dV*r;
   double Bz  = cons[BZZ]/dV;
   double B2 = Br*Br+Bp*Bp+Bz*Bz;

   double KE = .5*( Sr*vr + rho*vp_off*vp_off + Sz*vz );
   double rhoe = E-KE-.5*B2;
   double Pp = (gamma_law - 1.)*rhoe;

   if( Pp  < PRE_FLOOR*rho ) Pp = PRE_FLOOR*rho;
   if( isothermal ){
      double cs2 = get_cs2( x );
      Pp = cs2*rho/gamma_law;
   }

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr;
   prim[UPP] = vp/r;
   prim[UZZ] = vz;

   prim[BRR] = Br;
   prim[BPP] = Bp;
   prim[BZZ] = Bz;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      prim[q] = cons[q]/cons[DDD];
   }

}

void flux( const double * prim , double * flux , const double * x , const double * n, const double *xp, const double *xm ){
   UNUSED(xp);
   UNUSED(xm);

   double r = x[0];
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r;
   double vz  = prim[UZZ];

   double Br  = prim[BRR];
   double Bp  = prim[BPP];
   double Bz  = prim[BZZ];

   double vn = vr*n[0] + vp*n[1] + vz*n[2];
   double Bn = Br*n[0] + Bp*n[1] + Bz*n[2];
   double vB = vr*Br + vp*Bp + vz*Bz;

   double rhoe = Pp/(gamma_law-1.);
   double v2 = vr*vr + vp*vp + vz*vz;
   double B2 = Br*Br + Bp*Bp + Bz*Bz;

   flux[DDD] = rho*vn;
   flux[SRR] =     rho*vr*vn + (Pp+.5*B2)*n[0] - Br*Bn;
   flux[LLL] = r*( rho*vp*vn + (Pp+.5*B2)*n[1] - Bp*Bn );
   flux[SZZ] =     rho*vz*vn + (Pp+.5*B2)*n[2] - Bz*Bn;
   flux[TAU] = ( .5*rho*v2 + rhoe + Pp + B2 )*vn - vB*Bn;

   flux[BRR] =(Br*vn - vr*Bn) * bfield_scale_factor(r, 0);
   flux[BPP] =(Bp*vn - vp*Bn)/r;
   flux[BZZ] = Bz*vn - vz*Bn;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      flux[q] = prim[q]*flux[DDD];
   }
   
}

void source( const double * prim , double * cons , const double * xp , const double * xm , double dVdt ){
   
   UNUSED(xp);
   UNUSED(xm);

   double rp = xp[0];
   double rm = xm[0];
   double dphi = get_dp(xp[1],xm[1]);
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double r_1  = .5*(rp+rm);
   double r2_3 = (rp*rp + rp*rm + rm*rm)/3.;
   double vr  = prim[URR];
   double omega = prim[UPP];
   double r = get_centroid(rp, rm, 1);
   double z = get_centroid(xp[2], xm[2], 2);
   double x[3] = {r, 0.5*(xm[1]+xp[1]), z};

   double Br = prim[BRR];
   double Bp = prim[BPP];
   double Bz = prim[BZZ];

   double B2 = Br*Br+Bp*Bp+Bz*Bz;
 
   //double centrifugal = ( rho*omega*omega*r2_3 - Bp*Bp )/r_1*sin(.5*dphi)/(.5*dphi);
   double centrifugal;
   if(polar_sources)
   {
      centrifugal = ( rho*omega*omega*r2_3 )/r_1*sin(.5*dphi)/(.5*dphi);
      centrifugal -= Bp*Bp/r_1;
   }
   else
   {
       centrifugal = rho*omega*omega*r - Bp*Bp/r;
   }
   double press_bal   = (Pp+.5*B2)/r_1;

   cons[SRR] += dVdt*( centrifugal + press_bal );

   double om  = get_om( x );
   double om1 = get_om1( x );

   cons[TAU] += dVdt*rho*vr*( om*om*r2_3/r_1 - om1*(omega-om)*r2_3 );
 
   if( include_viscosity ){
      double nu = explicit_viscosity;
      cons[SRR] += -dVdt*nu*rho*vr/(r_1*r_1);
   }

}

void visc_flux(const double * prim, const double * gradr, const double * gradp,
               const double * gradz, double * flux,
               const double * x, const double * n)
{
   UNUSED(prim); 
   UNUSED(gradr); 
   UNUSED(gradp); 
   UNUSED(gradz); 
   UNUSED(flux); 
   UNUSED(x);
   UNUSED(n);
}

void visc_source(const double * prim, const double * gradr, const double *gradp,
                 const double * gradz, double * cons, const double *xp,
                 const double *xm, double dVdt)
{
   UNUSED(prim); 
   UNUSED(gradr); 
   UNUSED(gradp); 
   UNUSED(gradz); 
   UNUSED(cons); 
   UNUSED(xp);
   UNUSED(xm);
   UNUSED(dVdt);
}

void prim_to_E(const double *prim, double *E, const double *x)
{
    double r = x[0];

    double vr = prim[URR];
    double vp = prim[UPP]*r;
    double vz = prim[UZZ];
    double Br = prim[BRR];
    double Bp = prim[BPP];
    double Bz = prim[BZZ];

    E[0] = -(-vp*Bz+vz*Bp);
    E[1] = -vz*Br+vr*Bz;
    E[2] = -vr*Bp+vp*Br;
}

void flux_to_E( const double * Flux , const double * Ustr , const double * x , double * E1_riemann , double * B1_riemann , double * E2_riemann , double * B2_riemann , int dim ){

   double r = x[0];
   double irfac = 1.0/bfield_scale_factor(x[0], 0);

   if( dim==0 ){
       //PHI
      *E1_riemann =  Flux[BRR]*irfac;   // Ez 
      *B1_riemann =  Ustr[BRR]*irfac;   // Br
      *E2_riemann = -Flux[BZZ];         // Er 
      *B2_riemann =  Ustr[BZZ];         // Bz
   }else if( dim==1 ){
       //R
      *E1_riemann = -Flux[BPP]*r;       // Ez 
      *B1_riemann =  Ustr[BRR]*irfac;   // Br
      *E2_riemann =  Flux[BZZ];         // Ephi
   }else{
       //Z
      *E1_riemann =  Flux[BPP]*r;       // Er 
      *B1_riemann =  Ustr[BZZ];         // Bz
      *E2_riemann = -Flux[BRR]*irfac;   // Ephi
   }

}

void vel( const double * prim1 , const double * prim2 , double * Sl , double * Sr , double * Ss , const double * n , const double * x , double * Bpack ){
/*
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
*/

   double r = x[0];
   double L_Mins, L_Plus, L_Star;

   double P1   = prim1[PPP];
   double rho1 = prim1[RHO];

   double vr1  =   prim1[URR];
   double vp1  = r*prim1[UPP];
   double vz1  =   prim1[UZZ];

   double vn1  = vr1*n[0] + vp1*n[1] + vz1*n[2];

   double cs1  = sqrt( gamma_law*(P1/rho1) );

   double Br1 = prim1[BRR];
   double Bp1 = prim1[BPP];
   double Bz1 = prim1[BZZ];

   double Bn1 =  Br1*n[0] + Bp1*n[1] + Bz1*n[2];
   double B21 = (Br1*Br1  + Bp1*Bp1  + Bz1*Bz1);
   double b21 = B21/rho1;

   double FrL = vn1*Br1 - Bn1*vr1;
   double FpL = vn1*Bp1 - Bn1*vp1;
   double FzL = vn1*Bz1 - Bn1*vz1;

   double mrL = rho1*vr1;
   double mpL = rho1*vp1;
   double mzL = rho1*vz1;

   double FmrL = rho1*vr1*vn1 + (P1+.5*B21)*n[0] - Br1*Bn1;
   double FmpL = rho1*vp1*vn1 + (P1+.5*B21)*n[1] - Bp1*Bn1;
   double FmzL = rho1*vz1*vn1 + (P1+.5*B21)*n[2] - Bz1*Bn1;

   double cf21 = .5*( cs1*cs1 + b21 + sqrt(fabs(  (cs1*cs1+b21)*(cs1*cs1+b21) - 4.0*cs1*cs1*Bn1*Bn1/rho1 )) );

   L_Mins = vn1 - sqrt( cf21 );
   L_Plus = vn1 + sqrt( cf21 );

   double P2   = prim2[PPP];
   double rho2 = prim2[RHO];

   double vr2  =   prim2[URR];
   double vp2  = r*prim2[UPP];
   double vz2  =   prim2[UZZ];

   double vn2  = vr2*n[0] + vp2*n[1] + vz2*n[2];

   double cs2  = sqrt( gamma_law*(P2/rho2) );

   double Br2 = prim2[BRR];
   double Bp2 = prim2[BPP];
   double Bz2 = prim2[BZZ];

   double Bn2 =  Br2*n[0] + Bp2*n[1] + Bz2*n[2];
   double B22 = (Br2*Br2  + Bp2*Bp2  + Bz2*Bz2);
   double b22 = B22/rho2;

   double FrR = vn2*Br2 - Bn2*vr2;
   double FpR = vn2*Bp2 - Bn2*vp2;
   double FzR = vn2*Bz2 - Bn2*vz2;

   double mrR = rho2*vr2;
   double mpR = rho2*vp2;
   double mzR = rho2*vz2;

   double FmrR = rho2*vr2*vn2 + (P2+.5*B22)*n[0] - Br2*Bn2;
   double FmpR = rho2*vp2*vn2 + (P2+.5*B22)*n[1] - Bp2*Bn2;
   double FmzR = rho2*vz2*vn2 + (P2+.5*B22)*n[2] - Bz2*Bn2;

   double cf22 = .5*( cs2*cs2 + b22 + sqrt(fabs(  (cs2*cs2+b22)*(cs2*cs2+b22) - 4.0*cs2*cs2*Bn2*Bn2/rho2 )) );

   if( L_Mins > vn2 - sqrt( cf22 ) ) L_Mins = vn2 - sqrt( cf22 );
   if( L_Plus < vn2 + sqrt( cf22 ) ) L_Plus = vn2 + sqrt( cf22 );

   double aL = L_Plus;
   double aR = -L_Mins;

   double Br = ( aR*Br1 + aL*Br2 + FrL - FrR )/( aL + aR );
   double Bp = ( aR*Bp1 + aL*Bp2 + FpL - FpR )/( aL + aR );
   double Bz = ( aR*Bz1 + aL*Bz2 + FzL - FzR )/( aL + aR );
   double Bn = Br*n[0] + Bp*n[1] + Bz*n[2];

   double mr = ( aR*mrL + aL*mrR + FmrL - FmrR )/( aL + aR );
   double mp = ( aR*mpL + aL*mpR + FmpL - FmpR )/( aL + aR );
   double mz = ( aR*mzL + aL*mzR + FmzL - FmzR )/( aL + aR );

   double mnL = mrL*n[0]+mpL*n[1]+mzL*n[2];
   double mnR = mrR*n[0]+mpR*n[1]+mzR*n[2];
   double rho = ( aR*rho1 + aL*rho2 + mnL - mnR )/( aL + aR );

   L_Star = ( rho2*vn2*(L_Plus-vn2) - rho1*vn1*(L_Mins-vn1) + (P1+.5*B21-Bn1*Bn1) - (P2+.5*B22-Bn2*Bn2) )/( rho2*(L_Plus-vn2) - rho1*(L_Mins-vn1) );

   double vr = mr/rho;
   double vp = mp/rho;
   double vz = mz/rho;
   double vdotB = vr*Br + vp*Bp + vz*Bz;

   Bpack[0] = Bn;
   Bpack[1] = Br;
   Bpack[2] = Bp;
   Bpack[3] = Bz;
   Bpack[4] = vdotB;

   *Sl = L_Mins;
   *Sr = L_Plus;
   *Ss = L_Star;

}


double mindt(const double * prim , double w , const double * xp , const double * xm ){

   double r = get_centroid(xp[0], xm[0], 1);
   double Pp  = prim[PPP];
   double rho = prim[RHO];
   double vp  = (prim[UPP]-w)*r;
   double vr  = prim[URR];
   double vz  = prim[UZZ];
   double cs  = sqrt(gamma_law*Pp/rho);

   double Br  = prim[BRR];
   double Bp  = prim[BPP];
   double Bz  = prim[BZZ];
   double b2  = (Br*Br+Bp*Bp+Bz*Bz)/rho;

   double cf = sqrt(cs*cs + b2);

   double maxvr = cf + fabs(vr);
   double maxvp = cf + fabs(vp);
   double maxvz = cf + fabs(vz);

   double dtr = get_dL(xp,xm,1)/maxvr;
   double dtp = get_dL(xp,xm,0)/maxvp;
   double dtz = get_dL(xp,xm,2)/maxvz;
   
   double dt = dtr;
   if( dt > dtp ) dt = dtp;
   if( dt > dtz ) dt = dtz;

   return( dt );

}

double getReynolds( const double * prim , double w , const double * x , double dx ){

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

    if(dim == 0)
        return 1.0/x;
    else
        return 1.0;
}

double getCartInterpWeight(const double *x)
{
    UNUSED(x);
    return 0.0;
}
