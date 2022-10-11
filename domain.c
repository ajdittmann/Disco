
#include "paul.h"
#include "boundary.h"
#include "geometry.h"
#include "hydro.h"
#include "omega.h"
#include "analysis.h"
#include "planet.h"
#include "report.h"
#include "riemann.h"
#include "sink.h"
#include "output.h"

void setICparams( struct domain * );
void setHlldParams( struct domain * );
void setRotFrameParams( struct domain * );
void setMetricParams( struct domain * );
void setFrameParams(struct domain * );
void setNoiseParams( struct domain * );

int get_num_rzFaces( int , int , int );

void setupDomain( struct domain * theDomain ){

   srand(314159);
   rand();

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   theDomain->theCells = (struct cell **) malloc( Nr*Nz*sizeof(struct cell *) );
   int jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      theDomain->theCells[jk] = (struct cell *) malloc( Np[jk]*sizeof(struct cell) );
   }

   setGravParams( theDomain );
   setPlanetParams( theDomain );
   int Npl = theDomain->Npl;

   theDomain->thePlanets = NULL;
   theDomain->pl_gas_track = NULL;
   theDomain->pl_kin = NULL;
   theDomain->pl_RK_kin = NULL;
   theDomain->pl_aux = NULL;
   theDomain->pl_RK_aux = NULL;

   if(Npl > 0)
   {
       theDomain->thePlanets = (struct planet *) malloc(
                                    Npl * sizeof(struct planet) );

       theDomain->pl_gas_track = (double *) malloc(
                                    Npl * NUM_PL_INTEGRALS * sizeof(double));
       theDomain->pl_kin = (double *) malloc(
                                    Npl * NUM_PL_KIN * sizeof(double));
       theDomain->pl_RK_kin = (double *) malloc(
                                    Npl * NUM_PL_KIN * sizeof(double));
       theDomain->pl_aux = (double *) malloc(
                                    Npl * NUM_PL_AUX * sizeof(double));
       theDomain->pl_RK_aux = (double *) malloc(
                                    Npl * NUM_PL_AUX * sizeof(double));
   }
   setupPlanets(theDomain);

    //Setup Diagnostics & Snapshot

    setup_diagnostics(theDomain);
    setup_snapshot(theDomain);


    //Setup independent of node layout: pick the right rand()'s
    double Pmax = theDomain->phi_max;
    int j, k;

    //Discard everything from lower (global) z-layers.
    for(k=theDomain->N0z_glob; k<theDomain->N0z; k++)
        for(j=0; j<theDomain->Nr_glob; j++)
            rand();

    int i;

    for( k=0 ; k<Nz ; ++k )
    {
        //Discard randoms from inner (global) annuli
        for(j=theDomain->N0r_glob; j<theDomain->N0r; j++)
            rand();

        //DO the work
        for( j=0 ; j<Nr ; ++j )
        {
            jk = k*Nr + j;

            double p0 = Pmax*(double)rand()/(double)RAND_MAX;
            double dp = Pmax/(double)Np[jk];
            for( i=0 ; i<Np[jk] ; ++i )
            {
                double phi = p0+dp*(double)i;
                if( phi > Pmax ) phi -= Pmax;
                theDomain->theCells[jk][i].piph = phi;
                theDomain->theCells[jk][i].dphi = dp;
            }
        }
        //Discard randoms from outer (global) annuli
        for(j=theDomain->N0r+Nr; j<theDomain->N0r_glob + theDomain->Nr_glob;
                j++)
        {
            rand();
        }

    }

   theDomain->t       = theDomain->theParList.t_min;
   theDomain->t_init  = theDomain->theParList.t_min;
   theDomain->t_fin   = theDomain->theParList.t_max;

   theDomain->N_rpt = theDomain->theParList.NumRepts;
   theDomain->N_snp = theDomain->theParList.NumSnaps;
   theDomain->N_chk = theDomain->theParList.NumChecks;

   theDomain->count_steps = 0;
   theDomain->final_step  = 0;
   theDomain->check_plz   = 0;
   theDomain->startup = 1;

   theDomain->nrpt=-1;
   theDomain->nsnp=-1;
   theDomain->nchk=-1;

   theDomain->theFaces_1 = NULL;
   theDomain->theFaces_2 = NULL;
   theDomain->N_ftracks_r = get_num_rzFaces( Nr , Nz , 1 );
   theDomain->N_ftracks_z = get_num_rzFaces( Nr , Nz , 2 );

   theDomain->fIndex_r = (int *) malloc( (theDomain->N_ftracks_r+1)*sizeof(int) );
   theDomain->fIndex_z = (int *) malloc( (theDomain->N_ftracks_z+1)*sizeof(int) );

   setICparams( theDomain );
   setHydroParams( theDomain );
   setRiemannParams( theDomain );
   setHlldParams( theDomain );
   setOmegaParams( theDomain );
   setRotFrameParams( theDomain );
   setMetricParams( theDomain );
   setFrameParams( theDomain );
   setDiagParams( theDomain );
   setReportParams(theDomain);
   setNoiseParams( theDomain );
   setBCParams( theDomain );
   setSinkParams( theDomain );
}

void initial( double * , double * ); 
void restart( struct domain * ); 
void calc_dp( struct domain * );
void set_cell_xyz( struct domain * );
void set_wcell( struct domain * );
void adjust_gas( struct planet * , double * , double * , double );
void set_B_fields( struct domain * );
void subtract_omega( double * );
void addNoise(double *prim, double *x);
void exchangeData(struct domain *, int);

void setupCells( struct domain * theDomain ){

   int restart_flag = theDomain->theParList.restart_flag;
   if( restart_flag ) restart( theDomain );

   int noiseType = theDomain->theParList.noiseType;

   calc_dp( theDomain );
   set_cell_xyz(theDomain);

   int i,j,k;
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int NgRa = theDomain->NgRa;
   int NgRb = theDomain->NgRb;
   int NgZa = theDomain->NgZa;
   int NgZb = theDomain->NgZb;
   int Npl = theDomain->Npl;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int atmos = theDomain->theParList.include_atmos;

   //Null setup for all cells
   for(k=0; k<Nz; k++){
      for(j=0; j<Nr; j++){
         int jk = j+Nr*k;
         for(i=0; i<Np[jk]; i++){
            struct cell * c = &(theCells[jk][i]);
            c->wiph = 0.0; 
            c->real = 0;

            memset(c->gradr, 0, NUM_Q*sizeof(double));
            memset(c->gradp, 0, NUM_Q*sizeof(double));
            memset(c->gradz, 0, NUM_Q*sizeof(double));
         }
      }
   }

   //Setup real cells.
   for( k=NgZa ; k<Nz-NgZb ; ++k ){
      double z = get_centroid( z_kph[k], z_kph[k-1], 2);
      for( j=NgRa ; j<Nr-NgRb ; ++j ){
         double r = get_centroid( r_jph[j], r_jph[j-1], 1);
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = phip-c->dphi;
            c->wiph = 0.0; 
            double xp[3] = {r_jph[j  ],phip,z_kph[k  ]};
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            double dV = get_dV( xp , xm );
            double phi = c->piph-.5*c->dphi;
            double x[3] = {r, phi, z};

            if( !restart_flag )
            {
               initial( c->prim , x );
               subtract_omega( c->prim ); 
               if( atmos ){
                  int p;
                  for( p=0 ; p<Npl ; ++p ){
                     double gam = theDomain->theParList.Adiabatic_Index;
                     adjust_gas( theDomain->thePlanets+p , x , c->prim , gam );
                  }
               }
            }
            if(noiseType != 0)
                addNoise(c->prim, x);
            prim2cons( c->prim , c->cons , x , dV, xp, xm);
            c->real = 1;
         }    
      }    
   }

   if(!restart_flag && set_B_flag() && theDomain->theParList.CT)
   {
      // Communicate piph values to ghost zones.
      // TODO: WHY is this obly for CT??
      exchangeData(theDomain, 0);
      if( Nz > 1 )
         exchangeData(theDomain, 1);

      set_B_fields(theDomain);
   }

   for( k=NgZa ; k<Nz-NgZb ; ++k ){
      double z = get_centroid( z_kph[k], z_kph[k-1], 2);
      for( j=NgRa ; j<Nr-NgRb ; ++j ){
         double r = get_centroid( r_jph[j], r_jph[j-1], 1);
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = phip-c->dphi;
            double xp[3] = {r_jph[j  ],phip,z_kph[k  ]};
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            double dV = get_dV( xp , xm );
            double phi = c->piph-.5*c->dphi;
            double x[3] = {r, phi, z};
            cons2prim( c->cons , c->prim , x , dV, xp, xm);
         }
      }
   }

   set_wcell( theDomain );
}


/*
void clear_cell( struct cell * c ){
   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      c->prim[q]   = 0.0;
      c->cons[q]   = 0.0;
      c->RKcons[q] = 0.0;
      c->gradr[q]   = 0.0;
      c->gradp[q]  = 0.0;
      c->gradz[q]  = 0.0;
   }
   c->riph = 0.0;
   c->RKriph = 0.0;
   c->dr = 0.0;
   c->wiph = 0.0;
}
*/

void freeDomain( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      free( theDomain->theCells[jk] );
   }
   free( theDomain->theCells );
   free( theDomain->Np );
   theDomain->r_jph--;
   free( theDomain->r_jph );
   theDomain->z_kph--;
   free( theDomain->z_kph );

   if(theDomain->thePlanets != NULL)
      free( theDomain->thePlanets );
   if(theDomain->pl_gas_track != NULL)
      free( theDomain->pl_gas_track);
   if(theDomain->pl_kin != NULL)
      free( theDomain->pl_kin);
   if(theDomain->pl_RK_kin != NULL)
      free( theDomain->pl_RK_kin);
   if(theDomain->pl_aux != NULL)
      free( theDomain->pl_aux);
   if(theDomain->pl_RK_aux != NULL)
      free( theDomain->pl_RK_aux);


   free_diagnostics(&(theDomain->theTools));
   free_snapshot(&(theDomain->theSnap));

   free( theDomain->fIndex_r );
   free( theDomain->fIndex_z );

}

void check_dt( struct domain * theDomain , double * dt ){

   double t = theDomain->t;
   double tmax = theDomain->t_fin;
   int final=0;
   int check=0;
   if( t + *dt > tmax ){
      *dt = tmax-t;
      final=1;
   }

   if( theDomain->rank==0 ){
      FILE * abort = NULL;
      abort = fopen("abort","r");
      if( abort ){ final = 1; fclose(abort); }
      FILE * latest = NULL;
      latest = fopen("latest","r");
      if( latest ){ check = 1; fclose(latest); remove("latest");}
   }

#if USE_MPI
   MPI_Allreduce( MPI_IN_PLACE , &final , 1 , MPI_INT , MPI_SUM , theDomain->theComm );
   MPI_Allreduce( MPI_IN_PLACE , &check , 1 , MPI_INT , MPI_SUM , theDomain->theComm );
#endif
   if( final ) theDomain->final_step = 1;
   if( check ) theDomain->check_plz = 1;

}

void possiblyOutput( struct domain * theDomain , int override ){

   double t = theDomain->t;
   double t_min = theDomain->t_init;
   double t_fin = theDomain->t_fin;
   double Nrpt = theDomain->N_rpt;
   double Nsnp = theDomain->N_snp;
   double Nchk = theDomain->N_chk;
   int LogOut = theDomain->theParList.Out_LogTime;
   int restart_flag = theDomain->theParList.restart_flag;
   int n0;

   n0 = (int)( t*Nrpt/t_fin );
   if( LogOut ) n0 = (int)( Nrpt*log(t/t_min)/log(t_fin/t_min) );
   if( theDomain->nrpt < n0 || override ){
      theDomain->nrpt = n0;
      //longandshort( &theDomain , &L , &S , &iL , &iS , theDomain.theCells[0] , 0 , 0 );
      report( theDomain );
      //if( theDomain->rank==0 ) printf("t = %.5e\n",t);
   }
   n0 = (int)( t*Nchk/t_fin );
   if( LogOut ) n0 = (int)( Nchk*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nchk < n0 && Nchk>0) || override || theDomain->check_plz ){
      theDomain->nchk = n0;
      char filename[256];
      if( !override ){
         if( !theDomain->check_plz ){
            if(theDomain->rank==0)
                printf("Creating Checkpoint #%04d...\n",n0);
            if(theDomain->startup && restart_flag)
               sprintf(filename,"checkpoint_restart_%04d",n0);
            else
               sprintf(filename,"checkpoint_%04d",n0);
         }else{
            if(theDomain->rank==0) printf("Creating Requested Checkpoint...\n");
            sprintf(filename,"checkpoint_latest");
            theDomain->check_plz = 0;
         }
         output( theDomain , filename );
      }else{
         if(theDomain->rank==0) printf("Creating Final Checkpoint...\n");
         output( theDomain , "output" );
      }
   }
   n0 = (int)( t*Nsnp/t_fin );
   if( LogOut ) n0 = (int)( Nsnp*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nsnp < n0 && Nsnp>0) || override ){
      theDomain->nsnp = n0;
      char filename[256];
      if(!override)
          sprintf( filename , "snapshot_%06d" , n0 );
      else
          sprintf( filename , "snapshot_output" );
      snapshot( theDomain , filename );
   }

}


