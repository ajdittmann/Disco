#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void get_rpz(double *, double *);
void get_vec_from_rpz(double *, double *, double *);
void get_vec_from_xyz(double *, double *, double *);
void get_vec_contravariant(double *, double *, double *);

void initial( double * prim , double * x ){

    double rpz[3];
    get_rpz(x, rpz);
   double r   = rpz[0];
   double phi = rpz[1];
   double z   = rpz[2];

   double Rl = 0.45;
   double B0 = 1e-10;//1e-6;//0.0001;
   double Om = 1.0;
   double P0 = 1.0;

   double omega = Om;
   double Pp = P0 + .5*Om*Om*r*r;

   double xl = r*cos(phi);
   double yl = r*sin(phi);
   double zl = z;

   double rl = sqrt(xl*xl+zl*zl);
   double xx = M_PI*rl/Rl;

   double Bp = B0*pow(sin(xx),2.)*pow(cos(M_PI*yl/Rl),2.)*sqrt(2.*rl/Rl);
   if( rl > Rl || fabs(yl) > .5*Rl ) Bp = 0.0; 


   double Vrpz[3] = {0.0, r*omega, 0.0};
   double Bxyz[3] = {-Bp*zl/rl, 0.0, Bp*xl/rl};
   double V[3], B[3];
   get_vec_from_rpz(x, Vrpz, V);
   get_vec_contravariant(x, V, V);
   get_vec_from_xyz(x, Bxyz, B);

   prim[RHO] = 1.0; 
   prim[PPP] = Pp;
   prim[URR] = V[0]; 
   prim[UPP] = V[1];
   prim[UZZ] = V[2]; 

   prim[BRR] = B[0];
   prim[BPP] = B[1];
   prim[BZZ] = B[2]; 

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}
