
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "paul.h"
#include "geometry.h"


int getN0( int drank , int dsize , int dnum ){
   int N0 = (dnum*drank)/dsize;
   return(N0);
}


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
   if(!theDomain->theParList.NoBC_Rmin && Num_R > 1)
   {
       Nr_glob += Ng;
       N0r_glob -= Ng;
   }
   if(!theDomain->theParList.NoBC_Rmax && Num_R > 1)
       Nr_glob += Ng;
   if(!theDomain->theParList.NoBC_Zmin && Num_Z > 1)
   {
       Nz_glob += Ng;
       N0z_glob -= Ng;
   }
   if(!theDomain->theParList.NoBC_Zmax && Num_Z > 1)
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
      if( LogZoning == 0 ){
         theDomain->r_jph[j] = Rmin + x*(Rmax-Rmin);
      }else if( LogZoning == 1 ){
         theDomain->r_jph[j] = Rmin*pow(Rmax/Rmin,x);
      }else{
         theDomain->r_jph[j] = R0*(pow(Rmax/R0,x)-1) + Rmin + (R0-Rmin)*x;
      }
   }
   double dz = (Zmax-Zmin)/(double)Num_Z;
   double z0 = Zmin + (double)N0z*dz;
   for( k=-1 ; k<Nz ; ++k ){
      theDomain->z_kph[k] = z0 + ((double)k+1.)*dz;
   }

   double dr0;
   if(LogZoning == 0)
       dr0 = (Rmax-Rmin) / (double) Num_R;
   else if(LogZoning == 1)
       dr0 = Rmin * (pow(Rmax/Rmin, 1.0/Num_R) - 1.0);
   else
       dr0 = R0*(pow(Rmax/R0,1.0/Num_R)-1) + (R0-Rmin)/Num_R;
   theDomain->dr0 = dr0;
      

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


