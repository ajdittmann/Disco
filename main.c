
#include "paul.h"
#include "profiler.h"

int mpiSetup( struct domain * , int , char *[] );
void setupGrid( struct domain * );
void timestep( struct domain * , double );
void setupCells( struct domain * );
//void regrid( struct domain * );
void boundary_trans( struct domain * , int );
void exchangeData( struct domain * , int );
double getmindt( struct domain * );
void calc_prim( struct domain * );

void read_par_file( struct domain * );
 

void setupDomain( struct domain * );
void freeDomain( struct domain * );
void check_dt( struct domain * , double * );
void possiblyOutput( struct domain * , int );

void initializeReport(struct domain *);


void print_welcome();

int main( int argc , char * argv[] ){
 
#if USE_MPI
   MPI_Init(&argc,&argv);
#endif
   struct domain theDomain = {0};
   struct profiler prof;
   theDomain.prof = &prof;
   start_clock( &theDomain ); 
   read_par_file( &theDomain );
  
   if(theDomain.rank==0)
      print_welcome();

   int error = mpiSetup(&theDomain,argc,argv);
   if( error==1 ) return(0);

   if(theDomain.rank==0) remove("abort");

   setupGrid( &theDomain );   
   setupDomain( &theDomain );
 
   setupCells( &theDomain );
/*
   if( theDomain.theParList.Initial_Regrid && !(theDomain.theParList.restart_flag) ) regrid( &theDomain );
*/

   if( theDomain.Nr > 1 ){
      exchangeData( &theDomain , 0 );
      if( !theDomain.theParList.R_Periodic)
         boundary_trans( &theDomain , 1);
   }
   if( theDomain.Nz > 1 ){
      exchangeData( &theDomain , 1 );
      if( !theDomain.theParList.Z_Periodic)
         boundary_trans( &theDomain , 2);
   }

   initializeReport(&theDomain);

   while( !(theDomain.final_step) ){
      
      prof_tick(&prof, PROF_DT);
      double dt = getmindt( &theDomain );
      check_dt( &theDomain , &dt );
      prof_tock(&prof, PROF_DT);
      
      prof_tick(&prof, PROF_OUTPUT);
      possiblyOutput( &theDomain , 0 );
      prof_tock(&prof, PROF_OUTPUT);
     
      if(theDomain.rank == 0)
        printf("t: %.6le    dt: %.6le\n", theDomain.t, dt);
      prof_tick(&prof, PROF_TIMESTEP);
      timestep( &theDomain , dt );
      prof_tock(&prof, PROF_TIMESTEP);

      theDomain.startup = 0;
   }

   possiblyOutput( &theDomain , 1 );
   generate_log( &theDomain );
#if USE_MPI
   MPI_Barrier(theDomain.theComm);
#endif
   freeDomain( &theDomain );
#if USE_MPI
   MPI_Finalize();
#endif

   return(0);

}
