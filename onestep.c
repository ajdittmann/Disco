
#include "paul.h"
#include "hydro.h"
#include "profiler.h"
#include "analysis.h"
#include "planet.h"

void AMR( struct domain * ); 
void move_BCs( struct domain * , double );

void clean_pi( struct domain * );
void set_wcell( struct domain * );

void adjust_RK_cons( struct domain * , double );

void move_cells( struct domain * , double );
void set_cell_xyz(struct domain * theDomain);
void calc_dp( struct domain * );
void calc_prim( struct domain * );
void calc_cons( struct domain * );
void B_faces_to_cells( struct domain * , int );

void setup_faces( struct domain * , int );
void plm_phi( struct domain * );
void plm_trans( struct domain * , struct face * , int , int );
void phi_flux( struct domain * , double dt );
void trans_flux( struct domain * , double dt , int );
void add_source( struct domain * , double dt );

void avg_Efields( struct domain * );
void update_B_fluxes( struct domain * , double );
void subtract_advective_B_fluxes( struct domain * );
void check_flipped( struct domain * , int );
void flip_fluxes( struct domain * , int );

void boundary_trans( struct domain * , int );
void exchangeData( struct domain * , int );
void exchangePlanets(struct domain *theDomain);

void dump_grid(struct domain *, char filename[]);

//int get_num_rzFaces( int , int , int );

void checkNaNs(struct domain *theDomain, char label[])
{
    int Nz = theDomain->Nz;
    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;

    int i, j, k, q;
    int count_p = 0;
    int count_c = 0;
    for(k=0; k<Nz; k++)
        for(j=0; j<Nr; j++)
        {
            int jk = j + Nr*k;
            for(i=0; i<Np[jk]; i++)
            {
                struct cell *c = &(theDomain->theCells[jk][i]);

                //int flag = 0;
                for(q=0; q<NUM_Q; q++)
                {
                    if(c->prim[q] != c->prim[q])
                    {
                        count_p++;
                        //flag = 1;

                    }
                    if(c->cons[q] != c->cons[q])
                    {
                        count_c++;
                        //flag = 1;
                    }
                }
                //if(flag)
                //    printf("  NaN action at k = %d, j = %d, i = %d\n", k,j,i);
            }
        }
    if(count_p > 0)
        printf("NaNs in prim @ %s!\n", label);
    if(count_c > 0)
        printf("NaNs in cons @ %s!\n", label);
}


void onestep( struct domain * theDomain , double RK , double dt , int first_step , int last_step , double global_dt ){

   int Nz = theDomain->Nz;
   int Nr = theDomain->Nr;
   int bflag = set_B_flag();

   prof_tick(theDomain->prof, PROF_STEP_INIT);
   
   if( first_step ) set_wcell( theDomain );

   adjust_RK_cons( theDomain , RK );
   adjust_RK_diag( theDomain , RK );
   adjustPlanetsRKkin( theDomain , RK );
   adjustPlanetsRKaux( theDomain , RK );
   initializePlanetTracking(theDomain);
   
   prof_tock(theDomain->prof, PROF_STEP_INIT);

   //Reconstruction
   prof_tick(theDomain->prof, PROF_RECON);
   prof_tick(theDomain->prof, PROF_RECON_P);
   
   plm_phi( theDomain );
   
   prof_tock(theDomain->prof, PROF_RECON_P);

   if( Nr > 1 ){
      prof_tick(theDomain->prof, PROF_RECON_R);
      setup_faces( theDomain , 1 );
      int Nfr = theDomain->fIndex_r[theDomain->N_ftracks_r];
      plm_trans(theDomain, theDomain->theFaces_1, Nfr, 1);
      prof_tock(theDomain->prof, PROF_RECON_R);
   }
   if( Nz > 1 ){
      prof_tick(theDomain->prof, PROF_RECON_Z);
      setup_faces( theDomain , 2 );
      int Nfz = theDomain->fIndex_z[theDomain->N_ftracks_z];
      plm_trans(theDomain, theDomain->theFaces_2, Nfz, 2);
      prof_tock(theDomain->prof, PROF_RECON_Z);
   }
   prof_tock(theDomain->prof, PROF_RECON);

   //dump_grid(theDomain, "grid");

   //Flux
   prof_tick(theDomain->prof, PROF_FLUX);
   prof_tick(theDomain->prof, PROF_FLUX_P);

   phi_flux( theDomain , dt );
   prof_tock(theDomain->prof, PROF_FLUX_P);
   if( Nr > 1)
   {
      prof_tick(theDomain->prof, PROF_FLUX_R);
      trans_flux( theDomain , dt , 1 );
      prof_tock(theDomain->prof, PROF_FLUX_R);
   }
   if( Nz > 1 )
   {
      prof_tick(theDomain->prof, PROF_FLUX_Z);
      trans_flux( theDomain , dt , 2 );
      prof_tock(theDomain->prof, PROF_FLUX_Z);
   }
   prof_tock(theDomain->prof, PROF_FLUX);

   //CT update
   if( bflag && NUM_EDGES >= 4 ){
      prof_tick(theDomain->prof, PROF_CT);
      avg_Efields( theDomain );
      subtract_advective_B_fluxes( theDomain );
      update_B_fluxes( theDomain , dt );
      prof_tock(theDomain->prof, PROF_CT);
   }
  
   //Soucres
   prof_tick(theDomain->prof, PROF_SOURCE);
   add_source( theDomain , dt );
   prof_tock(theDomain->prof, PROF_SOURCE);

   prof_tick(theDomain->prof, PROF_MOVE);
   if( first_step ){
      move_cells( theDomain , dt );
      calc_dp( theDomain );
      set_cell_xyz(theDomain);
      if( bflag ){
         check_flipped( theDomain , 0 );
         flip_fluxes( theDomain , 0 );
         if( Nz>1 ){
            check_flipped( theDomain , 1 );
            flip_fluxes( theDomain , 1 );
         }
      }
   }
   prof_tock(theDomain->prof, PROF_MOVE);

   // Before we move the planets, update their internal diagnostics.
   // The gas_track integrals *must* be reduced (summed) over MPI ranks
   // now for LIVE planet motion.
   prof_tick(theDomain->prof, PROF_STEP_PL);
   if(!planet_motion_analytic())
      exchangePlanets(theDomain);

   updatePlanetsKinAux(theDomain, dt);

   if( !planet_motion_analytic()){
      movePlanetsLive(theDomain);
      setPlanetsXYZ(theDomain);
   }
   else if(first_step ){
      movePlanets( theDomain->thePlanets , theDomain->t , dt );
      setPlanetsXYZ(theDomain);
   }
   prof_tock(theDomain->prof, PROF_STEP_PL);

   prof_tick(theDomain->prof, PROF_CLEAN);
   clean_pi( theDomain );
   set_cell_xyz(theDomain);
   prof_tock(theDomain->prof, PROF_CLEAN);

   prof_tick(theDomain->prof, PROF_C2P);
   
   if( bflag && theDomain->theParList.CT ){
      B_faces_to_cells( theDomain , 1 );
   }
   
   calc_prim( theDomain ); //ORDERING??? AFTER?
   
   prof_tock(theDomain->prof, PROF_C2P);
   
   /*
   if( bflag && theDomain->theParList.CT ){
      B_faces_to_cells( theDomain , 0 );
   } 
   */

   //TODO: interaction with MHD? Hail Mary
   /*
   calc_cons(theDomain);
   */

   prof_tick(theDomain->prof, PROF_EXCHANGE);
   exchangeData( theDomain , 0 );
   prof_tock(theDomain->prof, PROF_EXCHANGE);

   if(! theDomain->theParList.R_Periodic)
   {
      prof_tick(theDomain->prof, PROF_BOUND);
      boundary_trans( theDomain , 1 );
      prof_tock(theDomain->prof, PROF_BOUND);
   }
   if( Nz > 1 ){
      prof_tick(theDomain->prof, PROF_EXCHANGE);
      exchangeData( theDomain , 1 );
      prof_tock(theDomain->prof, PROF_EXCHANGE);
      if(! theDomain->theParList.Z_Periodic)
      {
         prof_tick(theDomain->prof, PROF_BOUND);
         boundary_trans( theDomain , 2 );
         prof_tock(theDomain->prof, PROF_BOUND);
      }
   }
   
   //TODO: This was BEFORE BCs, but if wrecks cell pointers...
   //      Here, the BCs may not be satisfied if boundary zones are AMR'd...
   //TODO 2: AMR leading to STRANGE behaviour in 3d? Z boundaries being
   //           overwritten?  Needs a closer look.
   if( last_step ){
      //AMR( theDomain );
   }

   if( theDomain->theFaces_1 ) free( theDomain->theFaces_1 );
   if( theDomain->theFaces_2 ) free( theDomain->theFaces_2 );
}
