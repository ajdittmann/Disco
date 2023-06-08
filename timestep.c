#include "paul.h"
#include "analysis.h"
#include "planet.h"

void onestep( struct domain * , double , double , int , int);

void timestep( struct domain * theDomain , double dt ){
   
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int stepper = theDomain->theParList.Timestep;

   int i,jk;

   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         memcpy( c->RKcons , c->cons , NUM_Q*sizeof(double) );
#if NUM_FACES > 0
         memcpy( c->RK_Phi , c->Phi  , NUM_FACES*sizeof(double) );
#endif 
      }
   }

   copy_RK_diag(theDomain);
   copyPlanetsRK(theDomain);

   if(stepper == 1)
   {
      onestep( theDomain , 0.0 ,     dt , 1 , 1);
      theDomain->t += dt;
   }
   else if(stepper == 3)
   {
      onestep( theDomain ,   0.0,       dt, 1 , 0);
      theDomain->t += dt;
      onestep( theDomain ,  0.75,  0.25*dt, 0 , 0);
      theDomain->t -= dt*0.5;
      onestep( theDomain , 1./3., 2.*dt/3., 0 , 1);
      theDomain->t += dt*0.5;
   }
   else
   {
      onestep( theDomain , 0.0 ,     dt , 1 , 0);
      // Second RK2 timestep occurs at t^n+1
      theDomain->t += dt;   
      onestep( theDomain , 0.5 , 0.5*dt , 0 , 1);
   }

   add_diagnostics( theDomain , dt );

   theDomain->count_steps += 1;

}
