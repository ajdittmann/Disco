
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

double sStep(double x, double x0, double dx);

void setupGrid( struct domain * theDomain ){

   int Ng = NUM_G;
   theDomain->Ng = Ng;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   int Num_R = theDomain->theParList.Num_R;
   int Num_Z = theDomain->theParList.Num_Z;
   int LogZoning = theDomain->theParList.LogZoning;
   double aspect = theDomain->theParList.aspect;

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

   double R0 = theDomain->theParList.LogRadius;
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
         double xi;
         double drxi;
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
               drxi = a*pow(b, x);
            }
         }
         double dx = (1/factor)*dr/drxi;
         double a = 2*(1-factor)*dr/factor;
         double zval = sStep(0.0, xi, dx) + sStep(1.0, xi, dx)*x - sStep(x, xi, dx);
         theDomain->r_jph[j] += a*zval;
      }
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
/* Lambert W function. 
   Was ~/C/LambertW.c written K M Briggs Keith dot Briggs at bt dot com 97 May 21.  
   Revised KMB 97 Nov 20; 98 Feb 11, Nov 24, Dec 28; 99 Jan 13; 00 Feb 23; 01 Apr 09

   Computes Lambert W function, principal branch.
*/
   int i; 
   const double eps=4.0e-16, em1=0.3678794411714423215955237701614608; 
   double p,e,t,w;
   if (z<-em1 || isinf(z) || isnan(z)) { 
      fprintf(stderr,"LambertW: bad argument %g, exiting.\n",z); exit(1); 
   }
   if (0.0==z) return 0.0;
   if (z<-em1+1e-4) { // series near -em1 in sqrt(q)
      double q=z+em1,r=sqrt(q),q2=q*q,q3=q2*q;
      return 
       -1.0
       +2.331643981597124203363536062168*r
       -1.812187885639363490240191647568*q
       +1.936631114492359755363277457668*r*q
       -2.353551201881614516821543561516*q2
       +3.066858901050631912893148922704*r*q2
       -4.175335600258177138854984177460*q3
       +5.858023729874774148815053846119*r*q3
       -8.401032217523977370984161688514*q3*q;  // error approx 1e-16
   }
   /* initial approx for iteration... */
   if (z<1.0) { /* series near 0 */
      p=sqrt(2.0*(2.7182818284590452353602874713526625*z+1.0));
      w=-1.0+p*(1.0+p*(-0.333333333333333333333+p*0.152777777777777777777777)); 
   } else 
      w=log(z); /* asymptotic */
   if (z>3.0) w-=log(w); /* useful? */
   for (i=0; i<10; i++) { /* Halley iteration */
      e=exp(w); 
      t=w*e-z;
      p=w+1.0;
      t/=e*p-0.5*(p+1.0)*t/p; 
      w-=t;
      if (fabs(t)<eps*(1.0+fabs(w))) return w; /* rel-abs error */
   }
   /* should never get here */
   fprintf(stderr,"LambertW: No convergence at z=%g, exiting.\n",z); 
   exit(1);
}



double sStep(double x, double x0, double dx){
  return 0.5*(1.0 + tanh( (x-x0)/dx) );
}
