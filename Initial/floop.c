#include "../disco.h"

void get_xyz(double *x, double *xyz);
void get_rpz(double *x, double *rpz);
void get_vec_from_xyz(double *x, double *vxyz, double *v);
void get_vec_from_rpz(double *x, double *vrpz, double *v);
void get_vec_contravariant(double *x, double *v, double *vc);

static double *tp = NULL;

void setICparams( struct domain * theDomain ){
    tp = &(theDomain->t);
}

void initial( double * prim , double * x ){

   double xyz[3], rpz[3];
   get_xyz(x, xyz);
   get_rpz(x, rpz);
   double r = rpz[0];

   double Rl = 0.15;
   double B0 = 1e-4;
   double Om = 1.0;
   double P0 = 0.01;

   double x0 = 0.25;

   double omega = Om;
   double Pp = P0 + .5*Om*Om*r*r;

   double xl = xyz[0]-x0*cos(Om*(*tp));
   double yl = xyz[1]-x0*sin(Om*(*tp));

   double rl = sqrt(xl*xl+yl*yl);
   double xx = M_PI*rl/Rl;

   double Bp = B0*pow(sin(xx),2.)*sqrt(2.*rl/Rl);
   if( rl > Rl ) Bp = 0.0; 

   double dP = B0*B0*( - (rl/Rl)*pow(sin(xx),4.) - (1./16./M_PI)*( 12.*xx - 8.*sin(2.*xx) + sin(4.*xx) ) ); 
   if( rl > Rl ) dP = 0.0;

   double Vrpz[3] = {0.0, r*omega, 0.0};
   double Bxyz[3] = { -Bp*yl/rl, Bp*xl/rl, 0.0};
   double V[3], B[3];
   get_vec_from_xyz(x, Bxyz, B);
   get_vec_from_rpz(x, Vrpz, V);
   get_vec_contravariant(x, V, V);

   prim[RHO] = 1.0; 
   prim[PPP] = Pp+dP;
   prim[URR] = V[0]; 
   prim[UPP] = V[1];
   prim[UZZ] = V[2]; 

   prim[BRR] = B[0];
   prim[BPP] = B[1];
   prim[BZZ] = B[2]; 

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}
