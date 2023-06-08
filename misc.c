#include "paul.h"
#include "geometry.h"
#include "hydro.h"
#include "omega.h"
#include "planet.h"
#include "sink.h"
#include "riemann.h"

#include <string.h>

void clean_pi( struct domain * theDomain ){
   
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Npl = theDomain->Npl;
   double phi_max = theDomain->phi_max;

   int i,jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         double phi = c->piph;
         while( phi > phi_max ){ phi -= phi_max; }
         while( phi < 0.0     ){ phi += phi_max; }
         c->piph = phi; 
      }    
   }
   int p;
   for( p=0 ; p<Npl ; ++p ){
      struct planet * pl = theDomain->thePlanets + p;
      double phi = pl->phi;
      while( phi > phi_max ){ phi -= phi_max;}
      while( phi < 0.0     ){ phi += phi_max;}
      pl->phi = phi; 

      phi = theDomain->pl_kin[p*NUM_PL_KIN + PL_PHI];
      double RKphi = theDomain->pl_RK_kin[p*NUM_PL_KIN + PL_PHI];
      while( phi > phi_max ){ phi -= phi_max; RKphi -= phi_max; }
      while( phi < 0.0     ){ phi += phi_max; RKphi += phi_max; }
      theDomain->pl_kin[p*NUM_PL_KIN + PL_PHI] = phi;
      theDomain->pl_RK_kin[p*NUM_PL_KIN + PL_PHI] = RKphi;
   }

}


double getmindt( struct domain * theDomain ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int NgRa = theDomain->NgRa;
   int NgRb = theDomain->NgRb;
   int NgZa = theDomain->NgZa;
   int NgZb = theDomain->NgZb;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   double dt = theDomain->theParList.maxDT / theDomain->theParList.CFL;
   if(dt <= 0.0)
       dt = 1.0e100; //HUGE_VAL

   int i,j,k;
   for( k=NgZa ; k<Nz-NgZb ; ++k ){
      for( j=NgRa ; j<Nr-NgRb ; ++j ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            if(!(c->real))
                continue;
            double phip = c->piph;
            double phim = phip - c->dphi;
            double xp[3] = {r_jph[j  ] , phip , z_kph[k  ]};
            double xm[3] = {r_jph[j-1] , phim , z_kph[k-1]};
            int im = i-1;
            if( i==0 ) im = Np[jk]-1; 
            double wm = theCells[jk][im].wiph;
            double wp = c->wiph;
            double w = .5*(wm+wp);
            double dt_temp = mindt( c->prim , w , xp , xm );
            if( dt > dt_temp ) dt = dt_temp;
         }
      }
   }
   dt *= theDomain->theParList.CFL; 
#if USE_MPI
   MPI_Allreduce( MPI_IN_PLACE , &dt , 1 , MPI_DOUBLE , MPI_MIN , theDomain->theComm );
#endif

   return( dt );
}

void initial( double * , double * );
void restart( struct domain * );
/*
void clear_w( struct domain * theDomain ){
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int i,jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         theDomain->theCells[jk][i].wiph = 0.0;
      }
   }
}*/

void set_wcell( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int mesh_motion = theDomain->theParList.Mesh_Motion;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   int i,j,k;
   for( k=0 ; k<Nz ; ++k ){
      double z = get_centroid(z_kph[k], z_kph[k-1], 2);
      for( j=0 ; j<Nr ; ++j ){
         int jk = j+Nr*k;
         double r = get_centroid(r_jph[j], r_jph[j-1], 1);
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * cL = &(theCells[jk][i]);  
            double w = 0.0;
            if( mesh_motion ){
               double phip = cL->piph;
               double phim = phip-cL->dphi;
               double x[3] = {r, 0.5*(phim+phip), z};
               double wL = get_omega( cL->prim , x );

               int ip = (i+1)%Np[jk];
               struct cell * cR = &(theCells[jk][ip]);
               phip = cR->piph;
               phim = phip-cR->dphi;
               x[1] = 0.5*(phim+phip);
               double wR = get_omega( cR->prim , x );

               w = .5*(wL + wR); 
            }
            cL->wiph = w;
         }
      }    
   }
   if( mesh_motion == 3 ){
      for( k=0 ; k<Nz ; ++k ){
         for( j=0 ; j<Nr ; ++j ){
            int jk = j+Nr*k;
            double w = 0.0;
            for( i=0 ; i<Np[jk] ; ++i ){
               w += theCells[jk][i].wiph; 
            }    
            w /= (double)Np[jk];
            for( i=0 ; i<Np[jk] ; ++i ){
               theCells[jk][i].wiph = w; 
            }    
         }    
      } 
   }
   if( mesh_motion == 4 ){
      for( k=0 ; k<Nz ; ++k ){
         double z = get_centroid(z_kph[k], z_kph[k-1], 2);
         for( j=0 ; j<Nr ; ++j ){
            int jk = j+Nr*k;
            double r = get_centroid(r_jph[j], r_jph[j-1], 1);
            for( i=0 ; i<Np[jk] ; ++i ){
               double phip = theCells[jk][i].piph;
               double phim = phip-theCells[jk][i].dphi;
               double x[3] = {r, 0.5*(phim+phip), z};
               theCells[jk][i].wiph = mesh_om(x); 
            }    
         }    
      } 
   }
}

void clear_cell( struct cell * );

/*
void regrid( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   double jthresh = 0.1;
   int i,j,k;
   for( j=0 ; j<Nt ; ++j ){
   for( k=0 ; k<Np ; ++k ){
      int jk = j+Nt*k;
      for( i=2 ; i<Nr[jk]-1 ; ++i ){
         double jump = log( theCells[jk][i-1].prim[RHO]/theCells[jk][i+1].prim[RHO] );
         if( fabs(jump) > jthresh ){
            struct cell * c = theCells[jk]+i;
            double rp = c->riph;
            double rm = (c-1)->riph;
            int Nsplit = (int)(fabs(jump/jthresh));
            if(Nsplit>5) Nsplit=5;
            if(Nsplit>3) printf("r=%.2e jump=%.2e split = %d (No Cause For Alarm)\n",theCells[jk][i].riph,jump,Nsplit);
            int blocksize = (Nr[jk]-1) - i;
            Nr[jk] += Nsplit;
            theCells[jk] = (struct cell *) realloc( theCells[jk] , Nr[jk]*sizeof(struct cell) );
            memmove( theCells[jk]+i+1+Nsplit , theCells[jk]+i+1 , blocksize*sizeof(struct cell) );
            int l;
            for( l=0 ; l<Nsplit+1 ; ++l ){
               int m = l+i;
               struct cell * cnew = theCells[jk]+m;
               clear_cell( cnew );
               double rlp = rm + (rp-rm)*( (double)(l+1)/(double)(Nsplit+1) );
               double rlm = rm + (rp-rm)*( (double)(l  )/(double)(Nsplit+1) );
               double xp[3] = {rlp,t_jph[j]  ,p_kph[k]  };
               double xm[3] = {rlm,t_jph[j-1],p_kph[k-1]};
               double x[3];
               int d;
               for( d=0 ; d<3 ; ++d ) x[d] = .5*(xp[d]+xm[d]);
               initial( cnew->prim , x );
               cnew->riph = rlp;
               cnew->dr = rlp-rlm;
               double dV = get_dV(xp,xm);
               double rr = (2./3.)*(rlp*rlp*rlp-rlm*rlm*rlm)/(rlp*rlp-rlm*rlm);
               prim2cons( cnew->prim , cnew->cons , rr , dV );
            }
            i += Nsplit;
         }
      }
   }
   }
}
*/


void adjust_RK_cons( struct domain * theDomain , double RK ){
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;

   int i,jk,q;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         for( q=0 ; q<NUM_Q ; ++q ){
            c->cons[q] = (1.-RK)*c->cons[q] + RK*c->RKcons[q];
         }
#if NUM_FACES > 0
         for( q=0 ; q<NUM_FACES ; ++q ){
            c->Phi[q] = (1.-RK)*c->Phi[q] + RK*c->RK_Phi[q];
         }
#endif
      }
   }
}

void move_cells( struct domain * theDomain , double dt){
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int i,jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         c->piph += c->wiph*dt;
      }
   }
}

void set_cell_xyz(struct domain * theDomain)
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int * Np = theDomain->Np;
    int i,j,k;
   
    for(k=0; k<Nz; k++)
    {
        double zm = theDomain->z_kph[k-1];
        double zp = theDomain->z_kph[k];
        
        for(j=0; j<Nr; j++)
        {
            int jk = k*Nr + j;
            double rm = theDomain->r_jph[j-1];
            double rp = theDomain->r_jph[j];

            for( i=0 ; i<Np[jk] ; ++i )
            {
                struct cell *c = theDomain->theCells[jk] + i;
                double x[3];
                double xp[3] = {rp, c->piph, zp};
                double xm[3] = {rm, c->piph - c->dphi, zm};
                get_centroid_arr(xp, xm, x);
                get_xyz(x, c->xyz);
            }
        }
   }
}


void calc_dp( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;

   int i,jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         int im = i-1;
         if( i == 0 ) im = Np[jk]-1;
         double phim = theCells[jk][im].piph;
         double phip = theCells[jk][i ].piph;
         double dphi = get_dp(phip,phim);
         theCells[jk][i].dphi = dphi;
      }
   }
}

void calc_prim( struct domain * theDomain ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int NgRa = theDomain->NgRa;
   int NgRb = theDomain->NgRb;
   int NgZa = theDomain->NgZa;
   int NgZb = theDomain->NgZb;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   int i,j,k;
   for( k=NgZa ; k<Nz-NgZb ; ++k ){
      double zm = z_kph[k-1];
      double zp = z_kph[k];
      double z = get_centroid(zp, zm, 2);
      for( j=NgRa ; j<Nr-NgRb ; ++j ){
         int jk = j+Nr*k;
         double rm = r_jph[j-1];
         double rp = r_jph[j];
         double r = get_centroid(rp, rm, 1);
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = phip-c->dphi;
            double xp[3] = {rp, phip, zp};
            double xm[3] = {rm, phim, zm};
            double dV = get_dV( xp , xm );
            double x[3] = {r, 0.5*(phim+phip), z};
            cons2prim( c->cons , c->prim , x , dV, xp, xm );
         }
      }
   }
}

void calc_cons( struct domain * theDomain ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int NgRa = theDomain->NgRa;
   int NgRb = theDomain->NgRb;
   int NgZa = theDomain->NgZa;
   int NgZb = theDomain->NgZb;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   int i,j,k;
   for( k=NgZa ; k<Nz-NgZb ; ++k ){
      double zm = z_kph[k-1];
      double zp = z_kph[k];
      double z = get_centroid(zp, zm, 2);

      for( j=NgRa ; j<Nr-NgRb ; ++j ){
         double rm = r_jph[j-1];
         double rp = r_jph[j];
         double r = get_centroid(rp, rm, 1);

         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = phip-c->dphi;
            double xp[3] = {rp, phip, zp};
            double xm[3] = {rm, phim, zm};
            double dV = get_dV( xp , xm );
            double x[3] = {r, 0.5*(phim+phip), z};
            prim2cons( c->prim , c->cons , x , dV , xp, xm);
         }
      }
   }
}


void phi_flux( struct domain * theDomain , double dt ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int NgRa = theDomain->NgRa;
   int NgRb = theDomain->NgRb;
   int NgZa = theDomain->NgZa;
   int NgZb = theDomain->NgZb;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   int bflag = set_B_flag();
   int CT = theDomain->theParList.CT;
   int jmin, jmax, kmin, kmax;
   if(bflag && CT)
   {
       jmin = NgRa==0 ? 0  : NgRa-1;
       jmax = NgRb==0 ? Nr : Nr-NgRb+1;
       kmin = NgZa==0 ? 0  : NgZa-1;
       kmax = NgZb==0 ? Nz : Nz-NgZb+1;
   }
   else
   {
       jmin = NgRa;
       jmax = Nr-NgRb;
       kmin = NgZa;
       kmax = Nz-NgZb;
   }


   int i,j,k;
   for( k=kmin ; k<kmax ; ++k ){
      double zp = z_kph[k];
      double zm = z_kph[k-1];
      double z = get_centroid(zp, zm, 2);

      for( j=jmin ; j<jmax ; ++j ){
         double rp = r_jph[j];
         double rm = r_jph[j-1];
         double r = get_centroid(rp, rm, 1);
         
         int jk = j+Nr*k;
         struct cell * cp = theCells[jk];

         for( i=0 ; i<Np[jk] ; ++i ){
            int ip = i<Np[jk]-1 ? i+1 : 0;
            double phi = cp[i].piph;
            double xp[3] = {rp, phi, zp};
            double xm[3] = {rm, phi, zm};
            double dA = get_dA(xp,xm,0); 
            double x[3] = {r, phi, z};
            riemann_phi( &(cp[i]) , &(cp[ip]) , x , xp, xm, dA*dt );
         }
      }
   }

}

void buildfaces( struct domain * , int , int );

void trans_flux( struct domain * theDomain , double dt , int dim ){

    /*
     * Need to calculate fluxes for all interior faces.
     * For CT also need one level of fluxes outside the interior
     * ie. need the z-fluxes directly outside the r-boundary and vice-versa
     */ 

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int NgRa = theDomain->NgRa;
    int NgRb = theDomain->NgRb;
    int NgZa = theDomain->NgZa;
    int NgZb = theDomain->NgZb;
    
    int bflag = set_B_flag();
    int CT = theDomain->theParList.CT;

    int jmin, jmax, kmin, kmax, Nfr;

    int *fI;
    struct face * theFaces;

    double *fhydro_diag;
    double *fvisc_diag;

    if( dim==1 )
    {
        fI = theDomain->fIndex_r;
        theFaces = theDomain->theFaces_1;
        Nfr = Nr-1;
        jmin = NgRa==0 ? 0 : NgRa-1;
        jmax = NgRb==0 ? Nr-1 : Nr-NgRb;

        if(bflag && CT)
        {
            kmin = NgZa==0 ? 0 : NgZa-1;
            kmax = NgZb==0 ? Nz : Nz-NgZb+1;
        }
        else
        {
            kmin = NgZa;
            kmax = Nz-NgZb;
        }

        fhydro_diag = theDomain->theTools.F_r;
        fvisc_diag = theDomain->theTools.Fvisc_r;
    }
    else
    {
        fI = theDomain->fIndex_z;
        theFaces = theDomain->theFaces_2;
        Nfr = Nr;
        if(bflag && CT)
        {
            jmin = NgRa==0 ? 0  : NgRa-1;
            jmax = NgRb==0 ? Nr : Nr-NgRb+1;
        }
        else
        {
            jmin = NgRa;
            jmax = Nr-NgRb;
        }
        kmin = NgZa==0 ? 0  : NgZa-1;
        kmax = NgZb==0 ? Nz-1 : Nz-NgZb;
        
        fhydro_diag = theDomain->theTools.F_z;
        fvisc_diag = theDomain->theTools.Fvisc_z;
    }


    int j, k, q;
    for(k=kmin; k<kmax; k++)
    {
        double zm, zp;
        if(dim == 1)
            zm = theDomain->z_kph[k-1];
        else
            zm = theDomain->z_kph[k];
        zp = theDomain->z_kph[k];

        for(j=jmin; j<jmax; j++)
        {
            double rm, rp;
            if(dim == 1)
                rm = theDomain->r_jph[j];
            else
                rm = theDomain->r_jph[j-1];
            rp = theDomain->r_jph[j];

            int JK = j + Nfr*k;
            int f;
            for(f=fI[JK]; f<fI[JK+1]; f++)
            {
                double fdAdt_hydro[NUM_Q];
                double fdAdt_visc[NUM_Q];
                riemann_trans(theFaces + f, dt, dim, rp, rm, zp, zm,
                              fdAdt_hydro, fdAdt_visc);

                for(q=0; q<NUM_Q; q++)
                {
                    fhydro_diag[NUM_Q*JK+q] += fdAdt_hydro[q];
                    fvisc_diag[NUM_Q*JK+q] += fdAdt_visc[q];
                }
            }
        }
    }

   /*
   int f;
   for( f=0 ; f<Nf ; ++f ){
      riemann_trans( theFaces + f , dt , dim );
   }
   */
}


void setup_faces( struct domain * theDomain , int dim ){

   struct face ** theFaces;
   int NN;
   int * nn;
   if( dim==1 ){ 
      theFaces = &(theDomain->theFaces_1);
      nn = theDomain->fIndex_r;
      NN = theDomain->N_ftracks_r;
   }else{
      theFaces = &(theDomain->theFaces_2);
      nn = theDomain->fIndex_z;
      NN = theDomain->N_ftracks_z;
   }

   buildfaces( theDomain , dim , 0 );
   int Nf = nn[NN];
   *theFaces = (struct face *) malloc(((size_t) Nf) * sizeof(struct face));
   buildfaces( theDomain , dim , 1 );

}

void omega_src( double * , double * , double * , double * , double );

void add_source( struct domain * theDomain , double dt ){

   struct cell ** theCells = theDomain->theCells;
   struct planet * thePlanets = theDomain->thePlanets;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int NgRa = theDomain->NgRa;
   int NgRb = theDomain->NgRb;
   int NgZa = theDomain->NgZa;
   int NgZb = theDomain->NgZb;
   int * Np = theDomain->Np;
   int Npl = theDomain->Npl;

   int visc_flag = theDomain->theParList.visc_flag;

   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   int i,j,k,p;
   for( k=NgZa ; k<Nz-NgZb ; ++k ){
      for( j=NgRa ; j<Nr-NgRb ; ++j ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = c->piph - c->dphi;
            double xp[3] = {r_jph[j]  ,phip,z_kph[k]  };
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            double dV = get_dV(xp,xm);

            double sdVdt_hydro[NUM_Q] = {0};
            double sdVdt_grav[NUM_Q] = {0};
            double sdVdt_visc[NUM_Q] = {0};
            double sdVdt_sink[NUM_Q] = {0};
            double sdVdt_cool[NUM_Q] = {0};
            double sdVdt_damp[NUM_Q] = {0};
            
            source( c->prim , sdVdt_hydro, xp , xm , dV*dt  );
            
            for( p=0 ; p<Npl ; ++p ){
               planet_src( thePlanets+p, c->prim, sdVdt_grav, xp, xm, c->xyz,
                        dV, dt, theDomain->pl_gas_track + p*NUM_PL_INTEGRALS);
            }
            if(visc_flag)
                visc_source( c->prim, c->gradr, c->gradp, c->gradz, sdVdt_visc,
                            xp, xm, dV*dt);
            omega_src( c->prim , sdVdt_hydro , xp , xm , dV*dt );
            sink_src( c->prim , sdVdt_sink , xp , xm , c->xyz, dV, dt,
                      theDomain->pl_gas_track);
            cooling( c->prim , sdVdt_cool , xp , xm , c->xyz, dV, dt );
            damping( c->prim , sdVdt_damp , xp , xm , c->xyz, dV, dt );

            int q;
            for(q=0; q<NUM_Q; q++)
            {
                c->cons[q] += sdVdt_hydro[q] + sdVdt_grav[q] + sdVdt_visc[q]
                              + sdVdt_sink[q] + sdVdt_cool[q] + sdVdt_damp[q];
                int iq = NUM_Q*jk + q;
                theDomain->theTools.S[iq] += sdVdt_hydro[q];
                theDomain->theTools.Sgrav[iq] += sdVdt_grav[q];
                theDomain->theTools.Svisc[iq] += sdVdt_visc[q];
                theDomain->theTools.Ssink[iq] += sdVdt_sink[q];
                theDomain->theTools.Scool[iq] += sdVdt_cool[q];
                theDomain->theTools.Sdamp[iq] += sdVdt_damp[q];
            }
         }    
      }    
   }   
}

void longandshort( struct domain * theDomain , double * L , double * S , int * iL , int * iS , struct cell * sweep , int j , int k ){

   int Nr = theDomain->Nr;
   int jk = j + Nr*k;
   int Np = theDomain->Np[jk];
   double * r_jph = theDomain->r_jph;
   double dr = r_jph[j]-r_jph[j-1];
   double r  = r_jph[j];
   double Long = 0.0;
   double Short = 0.0;
   int iLong = 0;
   int iShort = 0;
   int i;
   for( i=0 ; i<Np ; ++i ){
      double dx = dr;
      double dy = r*sweep[i].dphi;
      double l = dy/dx;
      double s = dx/dy;
      if( Long  < l ){ Long  = l; iLong  = i; } 
      if( Short < s ){ Short = s; iShort = i; } 
   }
   *L  = Long;
   *iL = iLong;
   *S  = Short;
   *iS = iShort;

}

/*
void AMRsweep( struct domain * theDomain , struct cell ** swptr , int jk ){

   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int Nr = theDomain->Nr;
   int * Np = theDomain->Np;
   int j = jk%Nr;
   int k = jk/Nr;

   double MaxShort = theDomain->theParList.MaxShort;
   double MaxLong  = theDomain->theParList.MaxLong;

   struct cell * sweep = *swptr;
   double L,S;
   int iL=0;
   int iS=0;
   longandshort( theDomain , &L , &S , &iL , &iS , sweep , j , k );
   //printf("Long = %e, Short = %e\n",L,S);

   if( S>MaxShort ){
      int iSp = (iS+1)%Np[jk];

      //Possibly shift iS backwards by 1
      int iSm = iS-1;
      if( iSm == -1 ) iSm = Np[jk]-1;
      double dpL = sweep[iSm].dphi;
      double dpR = sweep[iSp].dphi;
      if( dpL < dpR ){
         --iS;
         --iSm;
         --iSp;
         if( iS  == -1 ) iS  = Np[jk]-1;
         if( iSm == -1 ) iSm = Np[jk]-1;
         if( iSp == -1 ) iSp = Np[jk]-1;
      }

      //Remove Zone at iS+1
      sweep[iS].dphi  += sweep[iSp].dphi;
      sweep[iS].piph   = sweep[iSp].piph;
      int q;
      for( q=0 ; q<NUM_Q ; ++q ){
         sweep[iS].cons[q]   += sweep[iSp].cons[q];
         sweep[iS].RKcons[q] += sweep[iSp].RKcons[q];
      }
      double phip = sweep[iS].piph;
      double phim = phip - sweep[iS].dphi;
      double xp[3] = {r_jph[j]  ,phip,z_kph[k]  };
      double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
      double r  = get_centroid( xp[0] , xm[0], 1 );
      double dV = get_dV( xp , xm );
      double x[3] = {r, 0.5*(phim+phip), 0.5*(z_kph[k-1]+z_kph[k])};
      cons2prim( sweep[iS].cons , sweep[iS].prim , x , dV, xp, xm );
      //Shift Memory
      int blocksize = Np[jk]-iSp-1;
      if( iSp != Np[jk]-1 ) memmove( sweep+iSp , sweep+iSp+1 , blocksize*sizeof(struct cell) );
      Np[jk] -= 1;
      *swptr = (struct cell *) realloc( sweep , Np[jk]*sizeof(struct cell) );
      sweep = *swptr;
      if( iS < iL ) iL--;

   }

   if( L>MaxLong ){
      Np[jk] += 1;
      *swptr = (struct cell *) realloc( sweep , Np[jk]*sizeof(struct cell) );
      sweep = *swptr;
      int blocksize = Np[jk]-iL-1;
      memmove( sweep+iL+1 , sweep+iL , blocksize*sizeof(struct cell) );

      double dphi = sweep[iL].dphi;
      double phip = sweep[iL].piph;
      double phim = phip - dphi;
      double phi0 = .5*(phip+phim);

      sweep[iL].piph   = phi0;
      sweep[iL].dphi   = .5*dphi;
      sweep[iL+1].dphi = .5*dphi;

      int q;
      for( q=0 ; q<NUM_Q ; ++q ){
         sweep[iL].cons[q]     *= .5;
         sweep[iL].RKcons[q]   *= .5;
         sweep[iL+1].cons[q]   *= .5;
         sweep[iL+1].RKcons[q] *= .5;
      }

      double xp[3] = {r_jph[j]  ,phi0,z_kph[k]  };
      double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
      double dV = get_dV( xp , xm );
      double r  = get_centroid( xp[0] , xm[0], 1 );
      double x[3] = {r, 0.5*(xp[1]+xm[1]), 0.5*(xp[2]+xm[2])};
      cons2prim( sweep[iL].cons , sweep[iL].prim , x , dV, xp, xm);

      xp[1] = phip;
      xm[1] = phi0;
      dV = get_dV( xp , xm );
      r  = get_centroid( xp[0] , xm[0], 1 );
      x[0] = r;
      x[1] = 0.5*(xp[1]+xm[1]);
      cons2prim( sweep[iL+1].cons , sweep[iL+1].prim , x , dV , xp, xm);

   }
}

void AMR( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      AMRsweep( theDomain , theCells+jk , jk );
   }
}
*/

void print_welcome()
{
    printf("\nDisco!\n\n");
    printf("Git Version: %s\n\n", GIT_VERSION);
    printf("*Compile-time Options*\n");
    printf("HYDRO: %s\n", HYDRO);
    printf("GEOMETRY: %s\n", GEOMETRY);
    printf("INITIAL: %s\n", INITIAL);
    printf("BOUNDARY: %s\n", BOUNDARY);
    printf("OUTPUT: %s\n", OUTPUT);
    printf("RESTART: %s\n", RESTART);
    printf("PLANETS: %s\n", PLANETS);
    printf("HLLD: %s\n", HLLD);
    printf("ANALYSIS: %s\n", ANALYSIS);
    printf("METRIC: %s\n", METRIC);
    printf("FRAME: %s\n", FRAME);
    printf("NUM_C: %d\n", NUM_C);
    printf("NUM_N: %d\n", NUM_N);
    printf("CT_MODE: %d\n", CT_MODE);
    printf("\n");
}


void dump_grid(struct domain *theDomain, char filename[])
{
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int NgRa = theDomain->NgRa;
   int NgRb = theDomain->NgRb;
   int NgZa = theDomain->NgZa;
   int NgZb = theDomain->NgZb;
   int * Np = theDomain->Np;

   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   int i,j,k, q;

   char fullName[256];
   sprintf(fullName, "%s.%05d.txt", filename, theDomain->count_steps);

   FILE *f = fopen(fullName, "w");

   for( k=NgZa ; k<Nz-NgZb ; ++k ){
      for( j=NgRa ; j<Nr-NgRb ; ++j ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = c->piph - c->dphi;
            double xp[3] = {r_jph[j]  ,phip,z_kph[k]  };
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            double x[3];
            get_centroid_arr(xp, xm, x);
            fprintf(f, "%d %d %d %.12le %.12le %.12le",
                    k, j, i, x[0], x[1], x[2]);
            for(q=0; q<NUM_Q; q++)
                fprintf(f, " %.12le", c->prim[q]);
            for(q=0; q<NUM_Q; q++)
                fprintf(f, " %.12le", c->cons[q]);
            for(q=0; q<NUM_Q; q++)
                fprintf(f, " %.12le", c->RKcons[q]);
            for(q=0; q<NUM_Q; q++)
                fprintf(f, " %.12le", c->gradr[q]);
            for(q=0; q<NUM_Q; q++)
                fprintf(f, " %.12le", c->gradp[q]);
            for(q=0; q<NUM_Q; q++)
                fprintf(f, " %.12le", c->gradz[q]);
            fprintf(f, "\n");
         }    
      }    
   }
   fclose(f);

}
