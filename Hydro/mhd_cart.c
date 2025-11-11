
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
   return(1);
}

double get_omega( const double * prim , const double * x ){
   return( prim[UPP] );
}


void prim2cons( const double * prim , double * cons , const double * x , double dV, const double *xp, const double *xm ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vx  = prim[URR];
   double vy  = prim[UPP];
   double vz  = prim[UZZ];
   double om  = get_om( x );
   double vy_off = vy - om;

   double Bx = prim[BRR];
   double By = prim[BPP];
   double Bz = prim[BZZ];

   double v2  = vx*vx + vy_off*vy_off + vz*vz;
   double B2  = Bx*Bx + By*By + Bz*Bz;

   double rhoe = Pp/(gamma_law - 1.);

   cons[DDD] = rho*dV;
   cons[TAU] = (.5*rho*v2 + rhoe + .5*B2 )*dV;
   cons[SRR] = rho*vx*dV;
   cons[LLL] = rho*vy*dV;
   cons[SZZ] = rho*vz*dV;

   cons[BRR] = Bx*dV;
   cons[BPP] = By*dV;
   cons[BZZ] = Bz*dV;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      cons[q] = prim[q]*cons[DDD];
   }
}

void getUstar( const double * prim , double * Ustar , const double * x , double Sk , double Ss , const double * n , const double * Bpack ){

   double Bsn = Bpack[0];
   double Bsx = Bpack[1];
   double Bsy = Bpack[2];
   double Bsz = Bpack[3];
   double vBs = Bpack[4];

   double rho = prim[RHO];
   double vx  = prim[URR];
   double vy  = prim[UPP];
   double vz  = prim[UZZ];
   double Pp  = prim[PPP];

   double Bx  = prim[BRR];
   double By  = prim[BPP];
   double Bz  = prim[BZZ];

   double v2 = vx*vx+vy*vy+vz*vz;
   double B2 = Bx*Bx+By*By+Bz*Bz;

   double vn = vx*n[0] + vy*n[1] + vz*n[2];
   double Bn = Bx*n[0] + By*n[1] + Bz*n[2];
   double vB = vx*Bx   + vy*By   + vz*Bz;

   double rhoe = Pp/(gamma_law-1.);

   double D  = rho; 
   double mx = rho*vx;
   double my = rho*vy;
   double mz = rho*vz;
   double E  = .5*rho*v2 + rhoe + .5*B2;

   double Bs2 = Bsx*Bsx+Bsy*Bsy+Bsz*Bsz;
   double Ps  = rho*( Sk - vn )*( Ss - vn ) + (Pp+.5*B2-Bn*Bn) - .5*Bs2 + Bsn*Bsn;

   double Dstar = ( Sk - vn )*D/( Sk - Ss );
   double Msx   = ( ( Sk - vn )*mx + ( Bx*Bn - Bsx*Bsn ) ) / ( Sk - Ss );
   double Msy   = ( ( Sk - vn )*my + ( By*Bn - Bsy*Bsn ) ) / ( Sk - Ss );
   double Msz   = ( ( Sk - vn )*mz + ( Bz*Bn - Bsz*Bsn ) ) / ( Sk - Ss );
   double Estar = ( ( Sk - vn )*E + (Ps+.5*Bs2)*Ss - (Pp+.5*B2)*vn - vBs*Bsn + vB*Bn ) / ( Sk - Ss );

   double Msn = Dstar*Ss;
   double mn  = Msx*n[0] + Msy*n[1] + Msz*n[2];

   Msx += n[0]*( Msn - mn );
   Msy += n[1]*( Msn - mn );
   Msz += n[2]*( Msn - mn );

   Ustar[DDD] = Dstar;
   Ustar[SRR] = Msx;
   Ustar[LLL] = Msy;
   Ustar[SZZ] = Msz;
   Ustar[TAU] = Estar;

   Ustar[BRR] = Bsx;
   Ustar[BPP] = Bsy;
   Ustar[BZZ] = Bsz;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      Ustar[q] = prim[q]*Ustar[DDD];
   }
}

void cons2prim( const double * cons , double * prim , const double * x , double dV, const double *xp, const double *xm  ){

   double rho = cons[DDD]/dV;
   if( rho < RHO_FLOOR )   rho = RHO_FLOOR;
   double Sx  = cons[SRR]/dV;
   double Sy  = cons[LLL]/dV;
   double Sz  = cons[SZZ]/dV;
   double E   = cons[TAU]/dV;
   double om  = get_om( x );
   
   double vx = Sx/rho;
   double vy = Sy/rho;
   double vy_off = vy - om;
   double vz = Sz/rho;

   double Bx  = cons[BRR]/dV;
   double By  = cons[BPP]/dV;
   double Bz  = cons[BZZ]/dV;
   double B2 = Bx*Bx+By*By+Bz*Bz;

   double KE = .5*( Sx*vx + rho*vy_off*vy_off + Sz*vz );
   double rhoe = E-KE-.5*B2;
   double Pp = (gamma_law - 1.)*rhoe;

   if( Pp  < PRE_FLOOR*rho ) Pp = PRE_FLOOR*rho;
   if( isothermal ){
      double cs2 = get_cs2( x );
      Pp = cs2*rho/gamma_law;
   }

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vx;
   prim[UPP] = vy;
   prim[UZZ] = vz;

   prim[BRR] = Bx;
   prim[BPP] = By;
   prim[BZZ] = Bz;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      prim[q] = cons[q]/cons[DDD];
   }
}

void flux( const double * prim , double * flux , const double * x , const double * n, const double *xp, const double *xm  ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vx  = prim[URR];
   double vy  = prim[UPP];
   double vz  = prim[UZZ];

   double Bx  = prim[BRR];
   double By  = prim[BPP];
   double Bz  = prim[BZZ];

   double vn = vx*n[0] + vy*n[1] + vz*n[2];
   double Bn = Bx*n[0] + By*n[1] + Bz*n[2];
   double vB = vx*Bx + vy*By + vz*Bz;

   double rhoe = Pp/(gamma_law-1.);
   double v2 = vx*vx + vy*vy + vz*vz;
   double B2 = Bx*Bx + By*By + Bz*Bz;

   flux[DDD] =   rho*vn;
   flux[SRR] =   rho*vx*vn + (Pp+.5*B2)*n[0] - Bx*Bn;
   flux[LLL] =   rho*vy*vn + (Pp+.5*B2)*n[1] - By*Bn;
   flux[SZZ] =   rho*vz*vn + (Pp+.5*B2)*n[2] - Bz*Bn;
   flux[TAU] = ( .5*rho*v2 + rhoe + Pp + B2 )*vn - vB*Bn;

   flux[BRR] = Bx*vn - vx*Bn;
   flux[BPP] = By*vn - vy*Bn;
   flux[BZZ] = Bz*vn - vz*Bn;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      flux[q] = prim[q]*flux[DDD];
   }
   
}

void source( const double * prim , double * cons , const double * xp , const double * xm , double dVdt ){
}

void visc_flux(const double * prim, const double * gradr, const double * gradp,
               const double * gradz, double * flux,
               const double * x, const double * n)
{

   int q;
   for(q=0; q<NUM_Q; q++)
      flux[q] = 0.0;
}

void visc_source(const double * prim, const double * gradr, const double *gradp,
                 const double * gradz, double * cons, const double *xp,
                 const double *xm, double dVdt){}

void prim_to_E(const double *prim, double *E, const double *x)
{
    double vx = prim[URR];
    double vy = prim[UPP];
    double vz = prim[UZZ];
    double Bx = prim[BRR];
    double By = prim[BPP];
    double Bz = prim[BZZ];

    E[0] = -vy*Bz+vz*By;
    E[1] = -vz*Bx+vx*Bz;
    E[2] = -vx*By+vy*Bx;
}

void flux_to_E( const double * Flux , const double * Ustr , const double * x , double * E1_riemann , double * B1_riemann , double * E2_riemann , double * B2_riemann , int dim ){

   if( dim==0 ){
       //Y
      *E1_riemann =  Flux[BRR]; // Ez 
      *B1_riemann =  Ustr[BRR]; // Bx
      *E2_riemann = -Flux[BZZ]; // Ex 
      *B2_riemann =  Ustr[BZZ]; // Bz
   }else if( dim==1 ){
       //X
      *E1_riemann = -Flux[BPP]; // Ez 
      *B1_riemann =  Ustr[BRR]; // Bx
      *E2_riemann =  Flux[BZZ]; // Ey
   }else{
       //Z
      *E1_riemann =  Flux[BPP]; // Ex 
      *B1_riemann =  Ustr[BZZ]; // Bz
      *E2_riemann = -Flux[BRR]; // Ey
   }
}

void vel( const double * prim1 , const double * prim2 , double * Sl , double * Sr , double * Ss , const double * n , const double * x , double * Bpack ){

   double L_Mins, L_Plus, L_Star;

   double P1   = prim1[PPP];
   double rho1 = prim1[RHO];

   double vx1  = prim1[URR];
   double vy1  = prim1[UPP];
   double vz1  = prim1[UZZ];

   double vn1  = vx1*n[0] + vy1*n[1] + vz1*n[2];

   double cs1  = sqrt( gamma_law*(P1/rho1) );

   double Bx1 = prim1[BRR];
   double By1 = prim1[BPP];
   double Bz1 = prim1[BZZ];

   double Bn1 =  Bx1*n[0] + By1*n[1] + Bz1*n[2];
   double B21 = (Bx1*Bx1  + By1*By1  + Bz1*Bz1);
   double b21 = B21/rho1;

   double FxL = vn1*Bx1 - Bn1*vx1;
   double FyL = vn1*By1 - Bn1*vy1;
   double FzL = vn1*Bz1 - Bn1*vz1;

   double mxL = rho1*vx1;
   double myL = rho1*vy1;
   double mzL = rho1*vz1;

   double FmxL = rho1*vx1*vn1 + (P1+.5*B21)*n[0] - Bx1*Bn1;
   double FmyL = rho1*vy1*vn1 + (P1+.5*B21)*n[1] - By1*Bn1;
   double FmzL = rho1*vz1*vn1 + (P1+.5*B21)*n[2] - Bz1*Bn1;

   double cf21 = .5*( cs1*cs1 + b21 + sqrt(fabs(  (cs1*cs1+b21)*(cs1*cs1+b21) - 4.0*cs1*cs1*Bn1*Bn1/rho1 )) );

   L_Mins = vn1 - sqrt( cf21 );
   L_Plus = vn1 + sqrt( cf21 );

   double P2   = prim2[PPP];
   double rho2 = prim2[RHO];

   double vx2  =   prim2[URR];
   double vy2  =   prim2[UPP];
   double vz2  =   prim2[UZZ];

   double vn2  = vx2*n[0] + vy2*n[1] + vz2*n[2];

   double cs2  = sqrt( gamma_law*(P2/rho2) );

   double Bx2 = prim2[BRR];
   double By2 = prim2[BPP];
   double Bz2 = prim2[BZZ];

   double Bn2 =  Bx2*n[0] + By2*n[1] + Bz2*n[2];
   double B22 = Bx2*Bx2  + By2*By2  + Bz2*Bz2;
   double b22 = B22/rho2;

   double FxR = vn2*Bx2 - Bn2*vx2;
   double FyR = vn2*By2 - Bn2*vy2;
   double FzR = vn2*Bz2 - Bn2*vz2;

   double mxR = rho2*vx2;
   double myR = rho2*vy2;
   double mzR = rho2*vz2;

   double FmxR = rho2*vx2*vn2 + (P2+.5*B22)*n[0] - Bx2*Bn2;
   double FmyR = rho2*vy2*vn2 + (P2+.5*B22)*n[1] - By2*Bn2;
   double FmzR = rho2*vz2*vn2 + (P2+.5*B22)*n[2] - Bz2*Bn2;

   double cf22 = .5*( cs2*cs2 + b22 + sqrt(fabs(  (cs2*cs2+b22)*(cs2*cs2+b22) - 4.0*cs2*cs2*Bn2*Bn2/rho2 )) );

   if( L_Mins > vn2 - sqrt( cf22 ) ) L_Mins = vn2 - sqrt( cf22 );
   if( L_Plus < vn2 + sqrt( cf22 ) ) L_Plus = vn2 + sqrt( cf22 );

   double aL = L_Plus;
   double aR = -L_Mins;

   double Bx = ( aR*Bx1 + aL*Bx2 + FxL - FxR )/( aL + aR );
   double By = ( aR*By1 + aL*By2 + FyL - FyR )/( aL + aR );
   double Bz = ( aR*Bz1 + aL*Bz2 + FzL - FzR )/( aL + aR );
   double Bn = Bx*n[0] + By*n[1] + Bz*n[2];

   double mx = ( aR*mxL + aL*mxR + FmxL - FmxR )/( aL + aR );
   double my = ( aR*myL + aL*myR + FmyL - FmyR )/( aL + aR );
   double mz = ( aR*mzL + aL*mzR + FmzL - FmzR )/( aL + aR );

   double mnL = mxL*n[0]+myL*n[1]+mzL*n[2];
   double mnR = mxR*n[0]+myR*n[1]+mzR*n[2];
   double rho = ( aR*rho1 + aL*rho2 + mnL - mnR )/( aL + aR );

   L_Star = ( rho2*vn2*(L_Plus-vn2) - rho1*vn1*(L_Mins-vn1) + (P1+.5*B21-Bn1*Bn1) - (P2+.5*B22-Bn2*Bn2) )/( rho2*(L_Plus-vn2) - rho1*(L_Mins-vn1) );

   double vx = mx/rho;
   double vy = my/rho;
   double vz = mz/rho;
   double vdotB = vx*Bx + vy*By + vz*Bz;

   Bpack[0] = Bn;
   Bpack[1] = Bx;
   Bpack[2] = By;
   Bpack[3] = Bz;
   Bpack[4] = vdotB;

   *Sl = L_Mins;
   *Sr = L_Plus;
   *Ss = L_Star;

}

double mindt(const double * prim , double w , const double * xp , const double * xm ){

   double Pp  = prim[PPP];
   double rho = prim[RHO];
   double vy  = (prim[UPP]-w);
   double vx  = prim[URR];
   double vz  = prim[UZZ];
   double cs  = sqrt(gamma_law*Pp/rho);

   double Bx  = prim[BRR];
   double By  = prim[BPP];
   double Bz  = prim[BZZ];
   double b2  = (Bx*Bx+By*By+Bz*Bz)/rho;

   double cf = sqrt(cs*cs + b2);

   double maxvx = cf + fabs(vx);
   double maxvy = cf + fabs(vy);
   double maxvz = cf + fabs(vz);

   double dtx = get_dL(xp,xm,1)/maxvx;
   double dty = get_dL(xp,xm,0)/maxvy;
   double dtz = get_dL(xp,xm,2)/maxvz;
   
   double dt = dtx;
   if( dt > dty ) dt = dty;
   if( dt > dtz ) dt = dtz;

   return( dt );

}

double getReynolds( const double * prim , double w , const double * x , double dx ){
   return 0;
}

void reflect_prims(double * prim, const double * x, int dim)
{
    //dim == 0: x, dim == 1: y, dim == 2: z
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
    // dim == 0: x, dim == 1: y, dim == 2: z

    return 1.0;
}

double getCartInterpWeight(const double *x)
{
    return 0.0;
}
