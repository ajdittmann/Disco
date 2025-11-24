
#include "../disco.h"
#include <string.h>

void initial( double * , double * );
double get_dV( double * , double * );
void cons2prim( double * , double * , double * , double );
void prim2cons( double * , double * , double * , double );

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

   double M = 1.0; //TODO: get from parfile
   double a = 1.0e1; //TODO: get from parfile
   double q = 1.0; //TODO: get from parfile

   double M1 = M / (1+q);
   double M2 = q * M / (1+q);
   double a1 = q * a / (1+q);
   double a2 = - a / (1+q);

   if( dim==1 && dim_rank[0] == dim_size[0]-1 ){
      int j;
      for( j=Nr-1 ; j>Nr-1-NgRb ; --j ){
         int i,k;
         for( k=0 ; k<Nz ; ++k ){
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i ){
               struct cell * c = &(theCells[jk][i]);
               double phi = c->piph - .5*c->dphi;
               double x[3] = { .5*(r_jph[j]+r_jph[j-1]) , phi , .5*(z_kph[k]+z_kph[k-1]) };
               initial( c->prim , x );
            }
         }
      }
   }

   if( dim==2 && dim_rank[1] == 0 ){
      int i,j,k;
      for( k=0 ; k<NgZa ; ++k ){
         for( j=0 ; j<Nr ; ++j ){
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i ){
               struct cell * c = &(theCells[jk][i]);
               double phi = c->piph - .5*c->dphi;
               double x[3] = { .5*(r_jph[j]+r_jph[j-1]) , phi , .5*(z_kph[k]+z_kph[k-1]) };
               initial( c->prim , x ); 
            }    
         }    
      } 
   }
   if( dim==2 && dim_rank[1] == dim_size[1]-1 ){ 
      int i,j,k;
      for( k=Nz-1 ; k>Nz-1-NgZb ; --k ){
         for( j=0 ; j<Nr ; ++j ){
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i ){
               struct cell * c = &(theCells[jk][i]);
               double phi = c->piph - .5*c->dphi;
               double x[3] = { .5*(r_jph[j]+r_jph[j-1]) , phi , .5*(z_kph[k]+z_kph[k-1]) };
               initial( c->prim , x );
            }
         }
      }
   }

   //Horizon jazz
   double rh1 = 2*M1;
   double rh2 = 2*M2;
   double cut1 = 0.8*rh1;
   double cut2 = 0.8*rh2;

   //Check if the cut-out region 1 is even on this process.
   if(r_jph[-1] < a1+cut1 && r_jph[Nr-1] > a1-cut1 
           && z_kph[-1] < cut1 && z_kph[Nz-1] > -cut1)
   {
       int i,j,k;
       for(k=0; k<Nz; k++)
       {
           double zp = z_kph[k];
           double zm = z_kph[k-1];
           double zc = fabs(zp) < fabs(zm) ? zp : zm;

           if(fabs(zc) > cut1 && Nz > 1)
               continue;

           for(j=0; j<Nr; j++)
           {
               double rp = r_jph[j];
               double rm = r_jph[j-1];
               double rc = fabs(rp-a1) < fabs(rm-a1) ? rp : rm;

               if(((rc-a1)*(rc-a1)+zc*zc > cut1*cut1 && Nz>1) || fabs(rc-a1) > cut1)
                   continue;

               int jk = j+Nr*k;
               for(i=0; i<Np[jk]; i++)
               {
                  struct cell * c = &(theCells[jk][i]);
                  double phi = c->piph - .5*c->dphi;
                  
                  double xc = rc*cos(phi);
                  double yc = rc*sin(phi);
                  double d2;
                  if(Nz > 1)
                    d2 = (xc-a1)*(xc-a1) + yc*yc + zc*zc;
                  else
                    d2 = (xc-a1)*(xc-a1) + yc*yc;

                  if(d2 < cut1*cut1)
                  {
                      double X[3] = {0.5*(rm+rp), phi, 0.5*(zm+zp)};
                      initial(c->prim, X);
                      c->real = 0;
                  }
               }
           }
       }
   }

   //Check if the cut-out region 2 is even on this process.
   if(r_jph[-1] < -a2+cut2 && r_jph[Nr-1] > -a2-cut2 
           && z_kph[-1] < cut2 && z_kph[Nz-1] > -cut2)
   {
       int i,j,k;
       for(k=0; k<Nz; k++)
       {
           double zp = z_kph[k];
           double zm = z_kph[k-1];
           double zc = fabs(zp) < fabs(zm) ? zp : zm;

           if(fabs(zc) > cut2 && Nz > 1)
               continue;

           for(j=0; j<Nr; j++)
           {
               double rp = r_jph[j];
               double rm = r_jph[j-1];
               double rc = fabs(rp+a2) < fabs(rm+a2) ? rp : rm;

               if(((rc+a2)*(rc+a2)+zc*zc > cut2*cut2 && Nz>1) || fabs(rc+a2) > cut2)
                   continue;

               int jk = j+Nr*k;
               for(i=0; i<Np[jk]; i++)
               {
                  struct cell * c = &(theCells[jk][i]);
                  double phi = c->piph - .5*c->dphi;
                  
                  double xc = rc*cos(phi);
                  double yc = rc*sin(phi);
                  double d2;
                  if(Nz > 1)
                    d2 = (xc-a2)*(xc-a2) + yc*yc + zc*zc;
                  else
                    d2 = (xc-a2)*(xc-a2) + yc*yc;

                  if(d2 < cut2*cut2)
                  {
                      double X[3] = {0.5*(rm+rp), phi, 0.5*(zm+zp)};
                      initial(c->prim, X);
                      c->real = 0;
                  }
               }
           }
       }
   }
}

