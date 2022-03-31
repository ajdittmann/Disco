#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "paul.h"
#include "geometry.h"


int getN0( int drank , int dsize , int dnum ){
   int N0 = (dnum*drank)/dsize;
   return(N0);
}


double LambertW(const double z);

void setupGrid( struct domain * theDomain ){

   int Ng = NUM_G;
   theDomain->Ng = Ng;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   int Num_R = theDomain->theParList.Num_R;
   int Num_Z = theDomain->theParList.Num_Z;
   int LogZoning = theDomain->theParList.LogZoning;
   double aspect = theDomain->theParList.aspect;

   double R0 = theDomain->theParList.LogRadius;
   double Rmin = theDomain->theParList.rmin;
   double Rmax = theDomain->theParList.rmax;
   double Zmin = theDomain->theParList.zmin;
   double Zmax = theDomain->theParList.zmax;
   double Pmax = theDomain->theParList.phimax;

   int Focus = theDomain->theParList.focusType;
   double focusPar1 = theDomain->theParList.focusPar1;
   double focusPar2 = theDomain->theParList.focusPar2;
   double focusPar3 = theDomain->theParList.focusPar3;
   double focusPar4 = theDomain->theParList.focusPar4;

   int N0r = getN0( dim_rank[0]   , dim_size[0] , Num_R );
   int N1r = getN0( dim_rank[0]+1 , dim_size[0] , Num_R );
   int NgRa = Ng;
   int NgRb = Ng;
   if( dim_rank[0] == 0 && theDomain->theParList.NoBC_Rmin)
       NgRa = 0;
   if( dim_rank[0] == dim_size[0]-1 && theDomain->theParList.NoBC_Rmax)
       NgRb = 0;
   N0r -= NgRa;
   N1r += NgRb;
   int Nr = N1r-N0r;

   int N0z = getN0( dim_rank[1]   , dim_size[1] , Num_Z );
   int N1z = getN0( dim_rank[1]+1 , dim_size[1] , Num_Z );
   int NgZa = Ng;
   int NgZb = Ng;
   if( Num_Z == 1 || (dim_rank[1] == 0 && theDomain->theParList.NoBC_Zmin))
       NgZa = 0;
   if( Num_Z == 1 || (dim_rank[1] == dim_size[1]-1
                        && theDomain->theParList.NoBC_Zmax))
       NgZb = 0;
   N0z -= NgZa;
   N1z += NgZb;
   int Nz = N1z-N0z;

   int Nr_glob = Num_R;
   int Nz_glob = Num_Z;
   int N0r_glob = 0;
   int N0z_glob = 0;
   if(!theDomain->theParList.NoBC_Rmin)
   {
       Nr_glob += Ng;
       N0r_glob -= Ng;
   }
   if(!theDomain->theParList.NoBC_Rmax)
       Nr_glob += Ng;
   if(!theDomain->theParList.NoBC_Zmin)
   {
       Nz_glob += Ng;
       N0z_glob -= Ng;
   }
   if(!theDomain->theParList.NoBC_Zmax)
       Nz_glob += Ng;

   theDomain->Nr = Nr;
   theDomain->Nz = Nz;
   theDomain->NgRa = NgRa;
   theDomain->NgRb = NgRb;
   theDomain->NgZa = NgZa;
   theDomain->NgZb = NgZb;
   theDomain->N0r = N0r;
   theDomain->N0z = N0z;
   theDomain->N0r_glob = N0r_glob;
   theDomain->N0z_glob = N0z_glob;
   theDomain->Nr_glob = Nr_glob;
   theDomain->Nz_glob = Nz_glob;
   printf("Rank = %d, Nr = %d, Nz = %d\n",theDomain->rank,Nr,Nz);

   theDomain->Np    = (int *)    malloc( Nr*Nz*sizeof(int) );
   theDomain->r_jph = (double *) malloc( (Nr+1)*sizeof(double) );
   theDomain->z_kph = (double *) malloc( (Nz+1)*sizeof(double) );

   ++(theDomain->r_jph);
   ++(theDomain->z_kph);

   int j,k;


   double xi, drxi, shift;
   double rshift = 0.0;
   double xs = 0.0;
   if(Focus > 0) {
      double r0 = R0;
      if (Focus==1) r0 = focusPar2;
      if (LogZoning == 0){
         xi = (r0-Rmin)/(Rmax-Rmin);
         drxi = Rmax-Rmin;
      }
      else if (LogZoning == 1){
         xi = log(r0/Rmin)/log(Rmax/Rmin);
         drxi = Rmin*log(Rmax/Rmin)*pow(Rmax/Rmin, xi);
      }else if (LogZoning == 2){
         double lC = (r0 - Rmin + R0)/(R0-Rmin);
         double lB = Rmax/R0;
         double lA = R0/(R0-Rmin);
         xi = (lC*log(lB) - LambertW(lA*pow(lB, lC)*log(lB)))/log(lB);
         drxi = R0-Rmin + R0*log(Rmax/R0)*pow(Rmax/R0, xi);
      }
      else {
         double x1 = 1.0/(1.0 - R0*log(R0/Rmax)/(R0-Rmin));
         double b = exp((1-Rmin/R0)/x1);
         double a = log(b)*Rmax/b;
         double c = a*pow(b,x1);
         if (r0 < R0) {
            xi = (r0-Rmin)/c;
            drxi = c;
         }else {
            xi = log(r0*log(b)/a)/log(b);
            drxi = a*pow(b, xi);
         }
      }
      if (Focus == 1) {
         double err = 10.0;
         double f;
         double shift = focusPar4;
         double factor = focusPar1;
         double dr = focusPar3;
         double dx = (1/factor)*dr/drxi;

         while (err > 1.e-10) {
            f = pow(cosh(rshift), -2.0) - pow(cosh(rshift - xi/dx), -2.0)*(1.0 - xi) - pow(cosh(rshift  + (1-xi)/dx), -2.0)*xi;
            f = (tanh(rshift) - tanh(rshift - xi/dx)*(1-xi) - tanh(rshift + (1-xi)/dx)*xi )/f;
            err = fabs(f);
            rshift -= f;
         }
         xs = xi  - shift*rshift*dx;
      }
   }

   for( j=-1 ; j<Nr ; ++j ){
      double x = (N0r + j + 1) / (double) Num_R;
      double value;
      if( LogZoning == 0 ){
         value = Rmin + x*(Rmax-Rmin);
      }else if( LogZoning == 1 ){
         value = Rmin*pow(Rmax/Rmin,x);
      }else if( LogZoning == 2){
         value = R0*(pow(Rmax/R0,x)-1) + Rmin + (R0-Rmin)*x;
      }else{
         double x1 = 1.0/(1.0 - R0*log(R0/Rmax)/(R0-Rmin));
         double b = exp((1-Rmin/R0)/x1);
         double a = log(b)*Rmax/b;
         double c = a*pow(b,x1);
         value = c*x + Rmin;
         if (x > x1) value = a*pow(b, x)/log(b);
      }
      theDomain->r_jph[j] = value;

      if(Focus == 1){
         double factor = focusPar1;
         double r0 = focusPar2;
         double dr = focusPar3;
         double shift = focusPar4;

         double dx = (1/factor)*dr/drxi;
         double a = 2*(1-factor)*dr/factor;

         double H = 0.5*(tanh((x-xs)/dx)+1);
         double H0 = 0.5*(tanh((-xs)/dx)+1);
         double H1 = 0.5*(tanh((1-xs)/dx)+1);
         theDomain->r_jph[j] += a*(H1*x + H0*(1-x) - H);
      }
   }
   for( j=0 ; j<Nr ; ++j ){
      if (theDomain->r_jph[j] <= theDomain->r_jph[j-1]) printf("YOUR GRID IS BAD at cell %d \n", j);
   }

   double dz = (Zmax-Zmin)/(double)Num_Z;
   double z0 = Zmin + (double)N0z*dz;
   for( k=-1 ; k<Nz ; ++k ){
      theDomain->z_kph[k] = z0 + ((double)k+1.)*dz;
   }

   theDomain->phi_max = theDomain->theParList.phimax;
   setGeometryParams( theDomain );

   for( k=0 ; k<Nz ; ++k ){
      double zp = theDomain->z_kph[k];
      double zm = theDomain->z_kph[k-1];

      for( j=0 ; j<Nr ; ++j ){
         int jk = j+Nr*k;
         double rp = theDomain->r_jph[j];
         double rm = theDomain->r_jph[j-1];

         double xp[3] = {rp, 0.0, zp};
         double xm[3] = {rm, 0.0, zm};
         double x[3];
         get_centroid_arr(xp, xm, x);

         double dr = get_dL(xp, xm, 1);
         double hp = get_scale_factor(xp, 0);

         double dp = aspect*dr/hp;
         int Np = (int)(Pmax/dp);
         if( Np<4 ) Np=4;
         theDomain->Np[jk] = Np;
      }
   }
}


double LambertW(const double z) {
   int i;
   const double eps=4.0e-16, em1=0.36787944117144232;
   double p,e,t,w;
   if (z < em1) {
      // Taylor series initial guess if within radius of convergence
      w = z - z*z + 1.5*pow(z,3) - 8*pow(z, 4)/3 + 125*pow(z,5)/24;
   } else {
      // otherwise use first few terms of asymptotic approximation
      double l1 = log(z);
      double l2 = log(log(z));
      w = l1 - l2 + l2/l1 + l2*(l2-2.0)*0.5/(l1*l1);
   }
   for (i=0; i<10; i++) {
      // Halley iteration, Corless et al. 1996, doi:10.1007/BF02124750
      e=exp(w);
      t=w*e-z;
      p=w+1.0;
      t = t/(e*p-(p+2)*t*0.5/p );
      w-=t;
      if (fabs(t)<eps*(1.0+fabs(w))) return w;
   }
   return w;
}

