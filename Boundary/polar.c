
#include "../disco.h"
#include <string.h>

void initial( double * , double * );
double get_dV( double * , double * );
void cons2prim( double * , double * , double * , double );
void prim2cons( double * , double * , double * , double );
double get_centroid( double , double , int);

void boundary_trans( struct domain * theDomain , int dim ){

   struct cell ** theCells = theDomain->theCells;

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int NgRb = theDomain->NgRb;
   int NgZa = theDomain->NgZa;
   int NgZb = theDomain->NgZb;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;

   if( dim==1 && dim_rank[0] == dim_size[0]-1 ){
      int j;
      for( j=Nr-1 ; j>Nr-1-NgRb ; --j ){
         double r = get_centroid( r_jph[j] , r_jph[j-1], 1);
         int i,k;
         for( k=0 ; k<Nz ; ++k ){
            double z = get_centroid( z_kph[k] , z_kph[k-1], 2);
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i ){
               struct cell * c = &(theCells[jk][i]);
               double phi = c->piph - .5*c->dphi;
               double x[3] = {r, phi, z};
               initial( c->prim , x );
            }
         }
      }
   }

   if( dim==2 && dim_rank[1] == 0 ){
      int i,j,k;
      for( k=0 ; k<NgZa ; ++k ){
         double z = get_centroid( z_kph[k] , z_kph[k-1], 2);
         for( j=0 ; j<Nr ; ++j ){
            int jk = j+Nr*k;
            double r = get_centroid( r_jph[j] , r_jph[j-1], 1 );
            for( i=0 ; i<Np[jk] ; ++i ){
               struct cell * c = &(theCells[jk][i]);
               double phi = c->piph - .5*c->dphi;
               double x[3] = {r, phi, z};
               initial( c->prim , x ); 
            }    
         }    
      } 
   }
   if( dim==2 && dim_rank[1] == dim_size[1]-1 ){ 
      int i,j,k;
      for( k=Nz-1 ; k>Nz-1-NgZb ; --k ){
         double z = get_centroid( z_kph[k] , z_kph[k-1], 2);
         for( j=0 ; j<Nr ; ++j ){
            double r = get_centroid( r_jph[j] , r_jph[j-1], 1 );
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i ){
               struct cell * c = &(theCells[jk][i]);
               double phi = c->piph - .5*c->dphi;
               double x[3] = { r , phi , z};
               initial( c->prim , x );
            }
         }
      }
   }

}

