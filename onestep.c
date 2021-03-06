
#include "paul.h"
#include "hydro.h"

void AMR( struct domain * ); 
void move_BCs( struct domain * , double );

void clean_pi( struct domain * );
void set_wcell( struct domain * );

void adjust_RK_cons( struct domain * , double );

void adjust_RK_planets_aux( struct domain * , double );
void adjust_RK_planets_kin( struct domain * , double );

void move_cells( struct domain * , double );
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

void movePlanets( struct planet * , double , double );
int planet_motion_analytic(void);

void boundary_trans( struct domain * , int );
void exchangeData( struct domain * , int );

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

   if( first_step ) set_wcell( theDomain );
   adjust_RK_cons( theDomain , RK );
   adjust_RK_planets_aux( theDomain , RK );

   //Reconstruction
   plm_phi( theDomain );
   if( Nr > 1 ){
      setup_faces( theDomain , 1 );
      int Nfr = theDomain->fIndex_r[theDomain->N_ftracks_r];
      plm_trans(theDomain, theDomain->theFaces_1, Nfr, 1);
   }
   if( Nz > 1 ){
      setup_faces( theDomain , 2 );
      int Nfz = theDomain->fIndex_z[theDomain->N_ftracks_z];
      plm_trans(theDomain, theDomain->theFaces_2, Nfz, 2);
   }

   //Flux
   phi_flux( theDomain , dt );
   if( Nr > 1)
      trans_flux( theDomain , dt , 1 );
   if( Nz > 1 )
      trans_flux( theDomain , dt , 2 );
   

   //CT update
   if( bflag && NUM_EDGES >= 4 ){
      avg_Efields( theDomain );
      subtract_advective_B_fluxes( theDomain );
      update_B_fluxes( theDomain , dt );
   }
  
   //Soucres
   add_source( theDomain , dt );

   if( first_step ){
      move_cells( theDomain , dt );
      if( bflag ){
         check_flipped( theDomain , 0 );
         flip_fluxes( theDomain , 0 );
         if( Nz>1 ){
            check_flipped( theDomain , 1 );
            flip_fluxes( theDomain , 1 );
         }
      }
   }
   

   if( !planet_motion_analytic() || first_step ){
      adjust_RK_planets_kin( theDomain , RK );
      movePlanets( theDomain->thePlanets , theDomain->t , dt );
   }
   clean_pi( theDomain );
   calc_dp( theDomain );

   
   if( bflag && theDomain->theParList.CT ){
      B_faces_to_cells( theDomain , 1 );
   }
   
   
   calc_prim( theDomain ); //ORDERING??? AFTER?
   
   /*
   if( bflag && theDomain->theParList.CT ){
      B_faces_to_cells( theDomain , 0 );
   } 
   */

   //TODO: interaction with MHD? Hail Mary
   /*
   calc_cons(theDomain);
   */

   exchangeData( theDomain , 0 );
   if(! theDomain->theParList.R_Periodic)
      boundary_trans( theDomain , 1 );
   if( Nz > 1 ){
      exchangeData( theDomain , 1 );
      if(! theDomain->theParList.Z_Periodic)
         boundary_trans( theDomain , 2 );
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
