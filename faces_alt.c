
#include "paul.h"
#include "geometry.h"

int get_num_rzFaces( int Nr , int Nz , int dim ){
   if( dim==1 ) return( (Nr-1)*Nz );
   else return( (Nz-1)*Nr );
}

void addFace( struct face * theFaces , int n , struct cell * cL , struct cell * cR , double dxL , double dxR , double * xp , double * xm , int dim , int LRtype ){
   double dp = get_dp(xp[1],xm[1]);
   double phic = get_dp(xp[1],.5*dp);
   double r = dim==1 ? 0.5*(xp[0]+xm[0]) : get_centroid(xp[0], xm[0], 1);
   double z = dim==2 ? 0.5*(xp[2]+xm[2]) : get_centroid(xp[2], xm[2], 2);

   theFaces[n].cm[0] = r;
   theFaces[n].cm[1] = phic;
   theFaces[n].cm[2] = z;
   theFaces[n].L   = cL;
   theFaces[n].R   = cR;
   theFaces[n].dxL = dxL;
   theFaces[n].dxR = dxR;
   theFaces[n].dphi= dp;
   theFaces[n].dA  = get_dA(xp,xm,dim);

   int dim_trans = 3-dim;//4-2*dim;//3-dim;
   theFaces[n].dl = get_dL( xp , xm , dim_trans );//xp[dim_trans] - xm[dim_trans];

   theFaces[n].E = 0.0;
   theFaces[n].B = 0.0;
   theFaces[n].LRtype = LRtype;
   theFaces[n].flip_flag = 0;

}

void buildfaces( struct domain * theDomain , int dim , int mode ){
  
   struct cell ** theCells = theDomain->theCells;
   struct face * theFaces;
   int * ntj;
   if( dim==1 ){
      theFaces = theDomain->theFaces_1;
      ntj = theDomain->fIndex_r;
   }else{
      theFaces = theDomain->theFaces_2;
      ntj = theDomain->fIndex_z;
   }

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   double Pmax = theDomain->phi_max;
   int j,k; 

   int I0[Nr*Nz];

   for( k=0 ; k<Nz ; ++k ){
      for( j=0 ; j<Nr ; ++j ){
         int jk = j+Nr*k;
         int found=0;
         int quad_prev=0;

         int i;
         for( i=0 ; i<Np[jk] && !found ; ++i ){
            struct cell * c = theCells[jk]+i;
            double convert = 2.*M_PI/Pmax;
            double sn = sin(c->piph*convert);
            double cs = cos(c->piph*convert);
            if( sn>0. && cs>0. && quad_prev ){
               quad_prev = 0;
               found = 1;
               I0[jk] = i; 
            }else if( sn<0. && cs>0. ){
               quad_prev=1;
            }
         }
         if( !found ) I0[jk]=0;
      }
   }

   int n=0;
   int Nrmax = Nr-1;
   if( dim==2 ) Nrmax = Nr;
   int Nzmax = Nz-1;
   if( dim==1 ) Nzmax = Nz;

   if( mode==0 ){
      for( k=0 ; k<Nzmax ; ++k ){
         for( j=0 ; j<Nrmax ; ++j ){
            int JK  = j+Nrmax*k;
            int jk  = j+Nr*k;
            int jkp = j+1 + Nr*k;
            if( dim==2 ) jkp = j+Nr*(k+1);

            ntj[JK] = n;
            n += Np[jk]+Np[jkp];
         }
      }
      ntj[Nrmax*Nzmax] = n;
   }else{
      for( k=0 ; k<Nzmax ; ++k ){
         for( j=0 ; j<Nrmax ; ++j ){
            int jp = j+1;
            int kp = k+1;

            int jk  = j  + Nr*k;
            int jkp = jp + Nr*k;
            if( dim==2 ) jkp = j + Nr*kp;

            double dxL,dxR;
            if( dim==1 ){
               double rm = get_centroid(r_jph[j], r_jph[j-1], 1);
               double rp = get_centroid(r_jph[jp], r_jph[j], 1);
               dxL = r_jph[j] - rm;
               dxR = rp - r_jph[j];
               //dxL = .5*(r_jph[j]  - r_jph[j-1]);
               //dxR = .5*(r_jph[jp] - r_jph[j]  );
            }else{
               double zm = get_centroid(z_kph[k], z_kph[k-1], 2);
               double zp = get_centroid(z_kph[kp], z_kph[k], 2);
               dxL = z_kph[k] - zm;
               dxR = zp - z_kph[k];
               //dxL = .5*(z_kph[k]  - z_kph[k-1]);
               //dxR = .5*(z_kph[kp] - z_kph[k]  );
            }
            double xp[3] = {r_jph[j],0.0,z_kph[k  ]};
            double xm[3] = {r_jph[j],0.0,z_kph[k-1]};
            if( dim==2 ){ xm[0] = r_jph[j-1] ; xm[2] = z_kph[k] ; }

            int i  = I0[jk];
            int ip = I0[jkp];

            int f;
            for( f=0 ; f<Np[jk]+Np[jkp] ; ++f ){
 
               struct cell * cL = &(theCells[jk ][i] );
               struct cell * cR = &(theCells[jkp][ip]);
 
               double dphi = cL->piph - cR->piph;
               while( dphi > .5*Pmax ) dphi -= Pmax;
               while( dphi <-.5*Pmax ) dphi += Pmax;
 
               int LR = 0;
               if( dphi > 0. ) LR = 1;

               dphi = cL->piph-cL->dphi - cR->piph+cR->dphi;
               while( dphi > .5*Pmax ) dphi -= Pmax;
               while( dphi <-.5*Pmax ) dphi += Pmax;

               int LR_back = 0;
               if( dphi < 0. ) LR_back = 1;

               if( LR == 0 ){
                  xp[1] = cL->piph;
               }else{
                  xp[1] = cR->piph;
               }
               if( LR_back == 0 ){
                  xm[1] = cL->piph-cL->dphi;
               }else{
                  xm[1] = cR->piph-cR->dphi;
               }
               addFace( theFaces , n , cL , cR , dxL , dxR , xp , xm , dim , LR );
               ++n;
            
               if( LR == 0 ){
                  ++i;
                  if( i == Np[jk] ) i=0;
               }else{
                  ++ip;
                  if( ip== Np[jkp] ) ip=0;
               }
            }
         }
      }
   }

}

