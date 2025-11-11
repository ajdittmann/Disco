
#include "../paul.h"
#include "../hydro.h"
#include "../geometry.h"
#include "../omega.h"

static double gamma_law = 0.0; 
static double RHO_FLOOR = 0.0; 
static double PRE_FLOOR = 0.0; 
static int isothermal = 0;
static int polar_sources_r = 0;
static int polar_sources_th = 0;

void setHydroParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
   isothermal = theDomain->theParList.isothermal_flag;
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
   if(theDomain->theParList.NoBC_Rmin == 1)
       polar_sources_r = 1;
   if(theDomain->theParList.NoBC_Zmin == 1
        || theDomain->theParList.NoBC_Zmax == 1)
       polar_sources_th = 1;
}

int set_B_flag(void){
   return(1);
}

double get_omega( const double * prim , const double * x ){
   return( prim[UPP] );
}


void prim2cons( const double * prim , double * cons , const double * x , double dV, const double *xp, const double *xm ){

   double r = x[0];
   double sinth = sin(x[2]);
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r*sinth;
   double vt  = prim[UZZ]*r;
   double om  = get_om( x );
   double vp_off = vp - om*r*sinth;

   double Br = prim[BRR];
   double Bp = prim[BPP];
   double Bt = prim[BZZ];

   double v2  = vr*vr + vp_off*vp_off + vt*vt;
   double B2  = Br*Br + Bp*Bp + Bt*Bt;

   double rhoe = Pp/(gamma_law - 1.);

   cons[DDD] = rho*dV;
   cons[TAU] = (.5*rho*v2 + rhoe + .5*B2 )*dV;
   cons[SRR] = rho*vr*dV;
   cons[LLL] = r*sinth*rho*vp*dV;
   cons[SZZ] = r*rho*vt*dV;

   cons[BRR] = Br*dV * bfield_scale_factor(r, 0);
   cons[BPP] = Bp*dV/(r*sinth);
   cons[BZZ] = Bt*dV/r;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      cons[q] = prim[q]*cons[DDD];
   }
}

void getUstar( const double * prim , double * Ustar , const double * x , double Sk , double Ss , const double * n , const double * Bpack ){

   double r = x[0];
   double sinth = sin(x[2]);

   double Bsn = Bpack[0];
   double Bsr = Bpack[1];
   double Bsp = Bpack[2];
   double Bst = Bpack[3];
   double vBs = Bpack[4];

   double rho = prim[RHO];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r*sinth;
   double vt  = prim[UZZ]*r;
   double Pp  = prim[PPP];

   double Br  = prim[BRR];
   double Bp  = prim[BPP];
   double Bt  = prim[BZZ];

   double v2 = vr*vr+vp*vp+vt*vt;
   double B2 = Br*Br+Bp*Bp+Bt*Bt;

   double vn = vr*n[0] + vp*n[1] + vt*n[2];
   double Bn = Br*n[0] + Bp*n[1] + Bt*n[2];
   double vB = vr*Br   + vp*Bp   + vt*Bt;

   double rhoe = Pp/(gamma_law-1.);

   double D  = rho; 
   double mr = rho*vr;
   double mp = rho*vp;
   double mt = rho*vt;
   double E  = .5*rho*v2 + rhoe + .5*B2;

   double Bs2 = Bsr*Bsr+Bsp*Bsp+Bst*Bst;
   double Ps  = rho*( Sk - vn )*( Ss - vn ) + (Pp+.5*B2-Bn*Bn) - .5*Bs2 + Bsn*Bsn;

   double Dstar = ( Sk - vn )*D/( Sk - Ss );
   double Msr   = ( ( Sk - vn )*mr + ( Br*Bn - Bsr*Bsn ) ) / ( Sk - Ss );
   double Msp   = ( ( Sk - vn )*mp + ( Bp*Bn - Bsp*Bsn ) ) / ( Sk - Ss );
   double Mst   = ( ( Sk - vn )*mt + ( Bt*Bn - Bst*Bsn ) ) / ( Sk - Ss );
   double Estar = ( ( Sk - vn )*E + (Ps+.5*Bs2)*Ss - (Pp+.5*B2)*vn - vBs*Bsn + vB*Bn ) / ( Sk - Ss );

   double Msn = Dstar*Ss;
   double mn  = Msr*n[0] + Msp*n[1] + Mst*n[2];

   Msr += n[0]*( Msn - mn );
   Msp += n[1]*( Msn - mn );
   Mst += n[2]*( Msn - mn );

   Ustar[DDD] = Dstar;
   Ustar[SRR] = Msr;
   Ustar[LLL] = r*sinth*Msp;
   Ustar[SZZ] = r*Mst;
   Ustar[TAU] = Estar;

   Ustar[BRR] = Bsr * bfield_scale_factor(r, 0);
   Ustar[BPP] = Bsp/(r*sinth);
   Ustar[BZZ] = Bst/r;

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

   double Br  = cons[BRR]/(dV * bfield_scale_factor(r, 0));
   double Bp  = cons[BPP]*r*sinth/dV;
   double Bt  = cons[BZZ]*r/dV;
   double B2 = Br*Br+Bp*Bp+Bt*Bt;

   double KE = .5*( Sr*vr + rho*vp_off*vp_off + St*vt );
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
   prim[UPP] = vp/(r*sinth);
   prim[UZZ] = vt/r;

   prim[BRR] = Br;
   prim[BPP] = Bp;
   prim[BZZ] = Bt;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      prim[q] = cons[q]/cons[DDD];
   }

}

void flux( const double * prim , double * flux , const double * x , const double * n, const double *xp, const double *xm ){

   double r = x[0];
   double sinth = sin(x[2]);
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r*sinth;
   double vt  = prim[UZZ]*r;

   double Br  = prim[BRR];
   double Bp  = prim[BPP];
   double Bt  = prim[BZZ];

   double vn = vr*n[0] + vp*n[1] + vt*n[2];
   double Bn = Br*n[0] + Bp*n[1] + Bt*n[2];
   double vB = vr*Br + vp*Bp + vt*Bt;

   double rhoe = Pp/(gamma_law-1.);
   double v2 = vr*vr + vp*vp + vt*vt;
   double B2 = Br*Br + Bp*Bp + Bt*Bt;

   flux[DDD] = rho*vn;
   flux[SRR] =           rho*vr*vn + (Pp+.5*B2)*n[0] - Br*Bn;
   flux[LLL] = r*sinth*( rho*vp*vn + (Pp+.5*B2)*n[1] - Bp*Bn );
   flux[SZZ] =       r*( rho*vt*vn + (Pp+.5*B2)*n[2] - Bt*Bn );
   flux[TAU] = ( .5*rho*v2 + rhoe + Pp + B2 )*vn - vB*Bn;

   flux[BRR] =(Br*vn - vr*Bn) * bfield_scale_factor(r, 0);
   flux[BPP] =(Bp*vn - vp*Bn)/(r*sinth);
   flux[BZZ] =(Bt*vn - vt*Bn)/r;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      flux[q] = prim[q]*flux[DDD];
   }
   
}

void source( const double * prim , double * cons , const double * xp , const double * xm , double dVdt ){
   
   double rho = prim[RHO];
   double Pp  = prim[PPP];

   double r = get_centroid(xp[0], xm[0], 1);
   double th = get_centroid(xp[2], xm[2], 2);
   double sinth = sin(th);
   double costh = cos(th);
   double up = prim[UPP];
   double ut = prim[UZZ];

   double Br = prim[BRR];
   double Bp = prim[BPP];
   double Bt = prim[BZZ];

   double B2 = Br*Br+Bp*Bp+Bt*Bt;
 
   //Polar_Sources are the result of integrating the centripetal source term
   //in a cartesian frame, assuming rho and omega are constant. This leads to
   //better behaviour at r=0.
   //
   //TODO: IMPLEMENT THIS
   //
   //The naive source term (polar_sources==0), on the other hand, can exactly
   //cancel with gravitational source terms.
   //
   double centrifugal_r, centrifugal_th;
   if(polar_sources_r)
      centrifugal_r = r*sinth*sinth*rho*up*up - Bp*Bp/r + r*rho*ut*ut-Bt*Bt/r;
   else
      centrifugal_r = r*sinth*sinth*rho*up*up - Bp*Bp/r + r*rho*ut*ut-Bt*Bt/r;
   if(polar_sources_th)
      centrifugal_th = r*r*sinth*costh*rho*up*up - Bp*Bp*costh/sinth;
   else
      centrifugal_th = r*r*sinth*costh*rho*up*up - Bp*Bp*costh/sinth;

   double r1_2 = 0.5*(xp[0]+xm[0]);
   double r2_3 = (xp[0]*xp[0] + xp[0]*xm[0] + xm[0]*xm[0]) / 3.0;
   double r3_4 = (3.0*xp[0]*r2_3 + xm[0]*xm[0]*xm[0]) / 4.0;
   double sinth1_2 = sin(0.5*(xp[2]+xm[2]));
   double costh1_2 = cos(0.5*(xp[2]+xm[2]));

   double Ptot = Pp + 0.5*B2;

   double press_bal_r, press_bal_th;
   press_bal_r = 2*Ptot*r1_2/r2_3;
   press_bal_th = Ptot*costh1_2/sinth1_2 * (r3_4*r1_2/(r2_3*r2_3));

   cons[SRR] += dVdt*( centrifugal_r + press_bal_r );
   cons[SZZ] += dVdt*( centrifugal_th + press_bal_th );


   //TODO: IMPLEMENT THIS
   //double om  = get_om( x );
   //double om1 = get_om1( x );

   //cons[TAU] += dVdt*rho*vr*( om*om*r2_3/r_1 - om1*(omega-om)*r2_3 );
}

void visc_flux(const double * prim, const double * gradr, const double * gradp,
               const double * gradz, double * flux,
               const double * x, const double * n){}

void visc_source(const double * prim, const double * gradr, const double *gradp,
                 const double * gradz, double * cons, const double *xp,
                 const double *xm, double dVdt){}

void prim_to_E(const double *prim, double *E, const double *x)
{
    double r = x[0];
    double sinth = sin(x[2]);

    double vr = prim[URR];
    double vp = prim[UPP]*r*sinth;
    double vt = prim[UZZ]*r;
    double Br = prim[BRR];
    double Bp = prim[BPP];
    double Bt = prim[BZZ];

    E[0] = -(-vp*Bt+vt*Bp);
    E[1] = -vt*Br+vr*Bt;
    E[2] = -vr*Bp+vp*Br;
}

void flux_to_E( const double * Flux , const double * Ustr , const double * x , double * E1_riemann , double * B1_riemann , double * E2_riemann , double * B2_riemann , int dim ){

   double r = x[0];
   double sinth = sin(x[2]);
   double irfac = 1.0/bfield_scale_factor(x[0], 0);

   if( dim==0 ){
       //PHI
      *E1_riemann =  Flux[BRR]*irfac;   // Et 
      *B1_riemann =  Ustr[BRR]*irfac;   // Br
      *E2_riemann = -Flux[BZZ]*r;       // Er 
      *B2_riemann =  Ustr[BZZ]*r;       // Bt
   }else if( dim==1 ){
       //R
      *E1_riemann = -Flux[BPP]*r*sinth; // Et
      *B1_riemann =  Ustr[BRR]*irfac;   // Br
      *E2_riemann =  Flux[BZZ]*r;       // Ephi
   }else{
       //Th
      *E1_riemann =  Flux[BPP]*r*sinth; // Er 
      *B1_riemann =  Ustr[BZZ]*r;       // Bt
      *E2_riemann = -Flux[BRR]*irfac;   // Ephi
   }

}

void vel( const double * prim1 , const double * prim2 , double * Sl , double * Sr , double * Ss , const double * n , const double * x , double * Bpack ){

   double r = x[0];
   double sinth = sin(x[2]);
   double L_Mins, L_Plus, L_Star;

   double P1   = prim1[PPP];
   double rho1 = prim1[RHO];

   double vr1  =   prim1[URR];
   double vp1  = r*sinth*prim1[UPP];
   double vt1  = r*prim1[UZZ];

   double vn1  = vr1*n[0] + vp1*n[1] + vt1*n[2];

   double cs1  = sqrt( gamma_law*(P1/rho1) );

   double Br1 = prim1[BRR];
   double Bp1 = prim1[BPP];
   double Bt1 = prim1[BZZ];

   double Bn1 =  Br1*n[0] + Bp1*n[1] + Bt1*n[2];
   double B21 = (Br1*Br1  + Bp1*Bp1  + Bt1*Bt1);
   double b21 = B21/rho1;

   double FrL = vn1*Br1 - Bn1*vr1;
   double FpL = vn1*Bp1 - Bn1*vp1;
   double FtL = vn1*Bt1 - Bn1*vt1;

   double mrL = rho1*vr1;
   double mpL = rho1*vp1;
   double mtL = rho1*vt1;

   double FmrL = rho1*vr1*vn1 + (P1+.5*B21)*n[0] - Br1*Bn1;
   double FmpL = rho1*vp1*vn1 + (P1+.5*B21)*n[1] - Bp1*Bn1;
   double FmtL = rho1*vt1*vn1 + (P1+.5*B21)*n[2] - Bt1*Bn1;

   double cf21 = .5*( cs1*cs1 + b21 + sqrt(fabs(  (cs1*cs1+b21)*(cs1*cs1+b21) - 4.0*cs1*cs1*Bn1*Bn1/rho1 )) );

   L_Mins = vn1 - sqrt( cf21 );
   L_Plus = vn1 + sqrt( cf21 );

   double P2   = prim2[PPP];
   double rho2 = prim2[RHO];

   double vr2  =         prim2[URR];
   double vp2  = r*sinth*prim2[UPP];
   double vt2  =       r*prim2[UZZ];

   double vn2  = vr2*n[0] + vp2*n[1] + vt2*n[2];

   double cs2  = sqrt( gamma_law*(P2/rho2) );

   double Br2 = prim2[BRR];
   double Bp2 = prim2[BPP];
   double Bt2 = prim2[BZZ];

   double Bn2 =  Br2*n[0] + Bp2*n[1] + Bt2*n[2];
   double B22 = (Br2*Br2  + Bp2*Bp2  + Bt2*Bt2);
   double b22 = B22/rho2;

   double FrR = vn2*Br2 - Bn2*vr2;
   double FpR = vn2*Bp2 - Bn2*vp2;
   double FtR = vn2*Bt2 - Bn2*vt2;

   double mrR = rho2*vr2;
   double mpR = rho2*vp2;
   double mtR = rho2*vt2;

   double FmrR = rho2*vr2*vn2 + (P2+.5*B22)*n[0] - Br2*Bn2;
   double FmpR = rho2*vp2*vn2 + (P2+.5*B22)*n[1] - Bp2*Bn2;
   double FmtR = rho2*vt2*vn2 + (P2+.5*B22)*n[2] - Bt2*Bn2;

   double cf22 = .5*( cs2*cs2 + b22 + sqrt(fabs(  (cs2*cs2+b22)*(cs2*cs2+b22) - 4.0*cs2*cs2*Bn2*Bn2/rho2 )) );

   if( L_Mins > vn2 - sqrt( cf22 ) ) L_Mins = vn2 - sqrt( cf22 );
   if( L_Plus < vn2 + sqrt( cf22 ) ) L_Plus = vn2 + sqrt( cf22 );

   double aL = L_Plus;
   double aR = -L_Mins;

   double Br = ( aR*Br1 + aL*Br2 + FrL - FrR )/( aL + aR );
   double Bp = ( aR*Bp1 + aL*Bp2 + FpL - FpR )/( aL + aR );
   double Bt = ( aR*Bt1 + aL*Bt2 + FtL - FtR )/( aL + aR );
   double Bn = Br*n[0] + Bp*n[1] + Bt*n[2];

   double mr = ( aR*mrL + aL*mrR + FmrL - FmrR )/( aL + aR );
   double mp = ( aR*mpL + aL*mpR + FmpL - FmpR )/( aL + aR );
   double mt = ( aR*mtL + aL*mtR + FmtL - FmtR )/( aL + aR );

   double mnL = mrL*n[0]+mpL*n[1]+mtL*n[2];
   double mnR = mrR*n[0]+mpR*n[1]+mtR*n[2];
   double rho = ( aR*rho1 + aL*rho2 + mnL - mnR )/( aL + aR );

   L_Star = ( rho2*vn2*(L_Plus-vn2) - rho1*vn1*(L_Mins-vn1) + (P1+.5*B21-Bn1*Bn1) - (P2+.5*B22-Bn2*Bn2) )/( rho2*(L_Plus-vn2) - rho1*(L_Mins-vn1) );

   double vr = mr/rho;
   double vp = mp/rho;
   double vt = mt/rho;
   double vdotB = vr*Br + vp*Bp + vt*Bt;

   Bpack[0] = Bn;
   Bpack[1] = Br;
   Bpack[2] = Bp;
   Bpack[3] = Bt;
   Bpack[4] = vdotB;

   *Sl = L_Mins;
   *Sr = L_Plus;
   *Ss = L_Star;

}


double mindt(const double * prim , double w , const double * xp , const double * xm ){

   double r  = get_centroid(xp[0], xm[0], 1);
   double th = get_centroid(xp[2], xm[2], 2);
   double sinth = sin(th);
   double Pp  = prim[PPP];
   double rho = prim[RHO];
   double vp  = (prim[UPP]-w)*r*sinth;
   double vr  = prim[URR];
   double vt  = prim[UZZ]*r;
   double cs  = sqrt(gamma_law*Pp/rho);

   double Br  = prim[BRR];
   double Bp  = prim[BPP];
   double Bt  = prim[BZZ];
   double b2  = (Br*Br+Bp*Bp+Bt*Bt)/rho;

   double cf = sqrt(cs*cs + b2);

   double maxvr = cf + fabs(vr);
   double maxvp = cf + fabs(vp);
   double maxvt = cf + fabs(vt);

   double dtr = get_dL(xp,xm,1)/maxvr;
   double dtp = get_dL(xp,xm,0)/maxvp;
   double dtt = get_dL(xp,xm,2)/maxvt;
   
   double dt = dtr;
   if( dt > dtp ) dt = dtp;
   if( dt > dtt ) dt = dtt;

   return( dt );

}

double getReynolds( const double * prim , double w , const double * x , double dx )
{
    return 0.0;
}

void reflect_prims(double * prim, const double * x, int dim)
{
    //dim == 0: r, dim == 1: p, dim == 2: t
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
    // dim == 0: r, dim == 1: p, dim == 2: t

    return 1.0;
}

double getCartInterpWeight(const double *x)
{
    return 0.0;
}
