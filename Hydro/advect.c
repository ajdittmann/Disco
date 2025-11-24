
#include "../disco.h"
#include "../hydro.h"
#include "../geometry.h"
#include "../omega.h"

static double RHO_FLOOR = 0.0; 
static int polar_sources = 0;
static int Cartesian_Interp = 0;

void setHydroParams(struct domain *theDomain)
{
    RHO_FLOOR = theDomain->theParList.Density_Floor;
    Cartesian_Interp = theDomain->theParList.Cartesian_Interp;
    if(theDomain->theParList.NoBC_Rmin == 1)
        polar_sources = 1;
}

int set_B_flag(void){
   return(0);
}

double get_omega( const double * prim , const double * x ){
   return( prim[UPP] );
}

void prim2cons(const double * prim, double * cons, const double * x,
               double dV, const double * xp, const double *xm){

   double rho = prim[RHO];

   cons[DDD] = rho*dV;
   cons[TAU] = 0.0;
   cons[SRR] = 0.0;
   cons[LLL] = 0.0;
   cons[SZZ] = 0.0;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      cons[q] = prim[q]*cons[DDD];
   }
}

void getUstar( const double * prim , double * Ustar , const double * x ,
              double Sk , double Ss , const double * n , const double * Bpack ){

   double r = x[0];
   double rho = prim[RHO];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r;
   double vz  = prim[UZZ];
   double vn = vr*n[0] + vp*n[1] + vz*n[2];


   double rhostar = rho*(Sk - vn)/(Sk - Ss);

   Ustar[DDD] = rhostar;
   Ustar[SRR] = 0.0;
   Ustar[LLL] = 0.0;
   Ustar[SZZ] = 0.0;
   Ustar[TAU] = 0.0;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      Ustar[q] = prim[q]*Ustar[DDD];
   }

}

void cons2prim( const double * cons , double * prim , const double * x ,
                double dV, const double * xp, const double *xm ){
   
   double rho = cons[DDD]/dV;
   if( rho < RHO_FLOOR )   rho = RHO_FLOOR;

   prim[RHO] = rho;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      prim[q] = cons[q]/cons[DDD];
   }

}

void flux( const double * prim , double * flux , const double * x ,
            const double * n, const double *xp, const double *xm ){
  
   double r = x[0];
   double rho = prim[RHO];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r;
   double vz  = prim[UZZ];

   double vn = vr*n[0] + vp*n[1] + vz*n[2];

   double scal_fac_n = 1.0;
   if(xp != NULL && xm != NULL && Cartesian_Interp)
   {
       double dphi = xp[1] - xm[1];
       double r_fac = n[0] > 0.0 ? 2*sin(0.5*dphi)/dphi : 1.0;
       scal_fac_n = r_fac * n[0] + n[1] + n[2];
   }

   flux[DDD] = rho*vn * scal_fac_n;
   flux[SRR] = 0.0;
   flux[LLL] = 0.0;
   flux[SZZ] = 0.0;
   flux[TAU] = 0.0;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      flux[q] = prim[q]*flux[DDD];
   }
   
}

void source( const double * prim , double * cons , const double * xp , const double * xm , double dVdt){
   
}

void visc_flux(const double * prim, const double * gradr, const double * gradp,
               const double * gradz, double * flux,
               const double * x, const double * n)
{
}

void visc_source(const double * prim, const double * gradr, const double *gradp,
                 const double * gradz, double * cons, const double *xp,
                 const double *xm, double dVdt)
{
}

void flux_to_E( const double * Flux , const double * Ustr , const double * x, 
                double * E1_riemann , double * B1_riemann , 
                double * E2_riemann , double * B2_riemann , int dim ){

   //Silence is Golden.

}

void vel( const double * prim1 , const double * prim2 , 
         double * Sl , double * Sr , double * Ss , 
         const double * n , const double * x , double * Bpack ){

   double r = x[0];
   double vn1  = prim1[URR]*n[0] + prim1[UPP]*n[1]*r + prim1[UZZ]*n[2];
   double vn2  = prim2[URR]*n[0] + prim2[UPP]*n[1]*r + prim2[UZZ]*n[2];

   *Sr =  vn1;
   *Sl =  vn1;

   if( *Sr < vn2 ) *Sr = vn2;
   if( *Sl > vn2 ) *Sl = vn2;

   *Ss = 0.5*(*Sr + *Sl);
}

double mindt(const double * prim , double w ,
             const double * xp , const double * xm ){

   double r = get_centroid(xp[0], xm[0], 1);
   double vp  = (prim[UPP]-w)*r;
   double vr  = prim[URR];
   double vz  = prim[UZZ];

   double maxvr = fabs(vr);
   double maxvp = fabs(vp);
   double maxvz = fabs(vz);

   double dtr = get_dL(xp,xm,1)/maxvr;
   double dtp = get_dL(xp,xm,0)/maxvp;
   double dtz = get_dL(xp,xm,2)/maxvz;
   
   double dt = dtr;
   if( dt > dtp ) dt = dtp;
   if( dt > dtz ) dt = dtz;

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
