#include "../paul.h"
#include <hdf5.h>
#include <string.h>
#include "../analysis.h"
#include "../boundary.h"
#include "../omega.h"
#include "../hydro.h"
#include "../planet.h"
#include "../report.h"

void getH5dims( char * file , char * group , char * dset , hsize_t * dims ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );
   hid_t h5spc = H5Dget_space( h5dst );

   H5Sget_simple_extent_dims( h5spc , dims , NULL);

   H5Sclose( h5spc );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void readSimple( char * file , char * group , char * dset , void * data , hid_t type ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   H5Dread( h5dst , type , H5S_ALL , H5S_ALL , H5P_DEFAULT , data );

   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void readPatch( char * file , char * group , char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   hsize_t mdims[dim];
   hsize_t fdims[dim];

   hsize_t fstart[dim];
   hsize_t fstride[dim];
   hsize_t fcount[dim];
   hsize_t fblock[dim];

   int d;
   for( d=0 ; d<dim ; ++d ){
      mdims[d] = (hsize_t) loc_size[d];
      fdims[d] = (hsize_t) glo_size[d];

      fstart[d]  = (hsize_t) start[d];
      fstride[d] = 1;
      fcount[d]  = (hsize_t) loc_size[d];
      fblock[d]  = 1;
   }
   hid_t mspace = H5Screate_simple(dim,mdims,NULL);
   hid_t fspace = H5Screate_simple(dim,fdims,NULL);

   H5Sselect_hyperslab( fspace , H5S_SELECT_SET , fstart , fstride , fcount , fblock );

   H5Dread( h5dst , type , mspace , fspace , H5P_DEFAULT , data );

   H5Sclose( mspace );
   H5Sclose( fspace );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}


void Doub2Cell( double * Q , struct cell * c ){
   int q;
   for( q=0 ; q<NUM_Q ; ++q ) c->prim[q] = Q[q];
#if NUM_FACES > 0
   for( q=0 ; q<NUM_FACES ; ++q ) c->Phi[q] = Q[NUM_Q+q];
#endif
   c->piph = Q[NUM_Q+NUM_FACES];
}

int getN0( int , int , int );
void freeDomain( struct domain * );

int get_num_rzFaces( int , int , int );

void setICparams( struct domain * );
void setRiemannParams( struct domain * );
void setGravParams( struct domain * );
void setHlldParams( struct domain * );
void setRotFrameParams( struct domain * );
void setMetricParams( struct domain * );
void setFrameParams(struct domain * );
void setDiagParams( struct domain * );
void setNoiseParams( struct domain * );
void setSinkParams( struct domain * );

void restart( struct domain * theDomain ){

   //This code has not been bug-tested in 3D.
   //The ordering of indices in particular hasn't been checked
   //i.e. should I use jk = j + Nt*k or jk = j*Np + k?

   freeDomain( theDomain );

   int Ng = theDomain->Ng;
   int rank = theDomain->rank;
   int size = theDomain->size;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;

   int i,j,k;

   char filename[256];
   strcpy(filename,"input.h5");
   char group1[256];
   strcpy(group1,"Grid");
   char group2[256];
   strcpy(group2,"Data");

   hsize_t dims[3];

   if( rank==0 ) printf("Restarting from file...\n");

   int NUM_R,NUM_Z,NUM_R_Tot,NUM_Z_Tot;
   double tstart;
   //Read the time from "T" and get Nt and Np from 
   //the dimensions of "Index".  Broadcast to all ranks.
#if USE_MPI
   if(rank==0){
#endif
      readSimple( filename , group1 ,"T", &tstart , H5T_NATIVE_DOUBLE );
      getH5dims( filename , group1 ,"Index", dims );
      NUM_Z_Tot = (int) dims[0];
      NUM_R_Tot = (int) dims[1];

#if USE_MPI
   }
#endif

#if USE_MPI
   MPI_Bcast( &NUM_R_Tot  , 1 , MPI_INT    , 0 , theDomain->theComm );
   MPI_Bcast( &NUM_Z_Tot  , 1 , MPI_INT    , 0 , theDomain->theComm );
   MPI_Bcast( &tstart , 1 , MPI_DOUBLE , 0 , theDomain->theComm );
#endif


   NUM_R = NUM_R_Tot;
   NUM_Z = NUM_Z_Tot;
   if(!theDomain->theParList.NoBC_Rmin)  NUM_R -= Ng;
   if(!theDomain->theParList.NoBC_Rmax)  NUM_R -= Ng;
   if(NUM_Z_Tot > 1 && !theDomain->theParList.NoBC_Zmin)  NUM_Z -= Ng;
   if(NUM_Z_Tot > 1 && !theDomain->theParList.NoBC_Zmax)  NUM_Z -= Ng;
   
   theDomain->theParList.Num_R = NUM_R;
   theDomain->theParList.Num_Z = NUM_Z;
   theDomain->t = tstart;
  
   //The following is very similar to equivalent code in gridsetup.c
   //Now you're just doing the process over because you're restarting
   //from file. 
   int N0r = getN0( dim_rank[0]   , dim_size[0] , NUM_R );
   int N1r = getN0( dim_rank[0]+1 , dim_size[0] , NUM_R );
   int NgRa = Ng;
   int NgRb = Ng;
   if( dim_rank[0] == 0 && theDomain->theParList.NoBC_Rmin)
       NgRa = 0;
   if( dim_rank[0] == dim_size[0]-1 && theDomain->theParList.NoBC_Rmax)
       NgRb = 0;
   N0r -= NgRa;
   N1r += NgRb;
   if(!theDomain->theParList.NoBC_Rmin) {
       N0r += Ng;
       N1r += Ng;
   }
   int Nr = N1r-N0r;

   int N0z = getN0( dim_rank[1]   , dim_size[1] , NUM_Z );
   int N1z = getN0( dim_rank[1]+1 , dim_size[1] , NUM_Z );
   int NgZa = Ng;
   int NgZb = Ng;
   if( NUM_Z == 1 || (dim_rank[1] == 0 && theDomain->theParList.NoBC_Zmin))
       NgZa = 0;
   if( NUM_Z == 1 || (dim_rank[1] == dim_size[1]-1
                        && theDomain->theParList.NoBC_Zmax))
       NgZb = 0;
   N0z -= NgZa;
   N1z += NgZb;
   if(NUM_Z > 1 && !theDomain->theParList.NoBC_Zmin) {
       N0z += Ng;
       N1z += Ng;
   }
   int Nz = N1z-N0z;
   
   theDomain->Nr = Nr;
   theDomain->Nz = Nz;
   theDomain->NgRa = NgRa;
   theDomain->NgRb = NgRb;
   theDomain->NgZa = NgZa;
   theDomain->NgZb = NgZb;

   theDomain->Np    = (int *)    malloc( ((size_t) (Nr*Nz)) * sizeof(int) );
   theDomain->r_jph = (double *) malloc( ((size_t) (Nr+1)) * sizeof(double) );
   theDomain->z_kph = (double *) malloc( ((size_t) (Nz+1)) * sizeof(double) );

   theDomain->theCells = (struct cell **) malloc( ((size_t) (Nr*Nz))
                                                    * sizeof( struct cell * ) );

   struct cell ** theCells = theDomain->theCells;

   //The following must happen in serial because different processors
   //will try to read from the same file.  In principle we can write
   //this using parallel HDF5, but I'd rather cross that bridge 
   //when I come to it.
   int nrk;
   int Nq=0;
   for( nrk=0 ; nrk<size ; ++nrk ){
   if( rank==nrk ){

      //Read the R values of the grid...
      int start1[1] = {N0r};
      int loc_size1[1] = {Nr+1};
      int glo_size1[1] = {NUM_R_Tot+1};
      double r_jph[Nr+1];
      readPatch( filename , group1 ,"r_jph", r_jph , H5T_NATIVE_DOUBLE , 1 , start1 , loc_size1 , glo_size1 ); 
      memcpy( theDomain->r_jph , r_jph , ((size_t) (Nr+1)) * sizeof(double) );
 
      //Read the Z values of the grid...
      start1[0]    = N0z;
      loc_size1[0] = Nz+1;
      glo_size1[0] = NUM_Z_Tot+1;
      double z_kph[Nz+1];
      readPatch( filename , group1 ,"z_kph", z_kph , H5T_NATIVE_DOUBLE , 1 , start1 , loc_size1 , glo_size1 );
      memcpy( theDomain->z_kph , z_kph , ((size_t) (Nz+1)) * sizeof(double) );

      ++(theDomain->r_jph);
      ++(theDomain->z_kph);

      //Read the indexing information so you know how to read in
      //The radial tracks of data which are coming up...
      int start2[2]   = {N0z,N0r};
      int loc_size2[2] = {Nz,Nr};
      int glo_size2[2] = {NUM_Z_Tot,NUM_R_Tot};
      int Np[Nr*Nz];
      int Index[Nr*Nz];

      printf("%d %d\n", NUM_Z, NUM_R);
      printf("%d %d\n", NUM_Z_Tot, NUM_R_Tot);
      printf("%d %d\n", Nz, Nr);
      printf("%d %d\n", N0z, N0r);

      readPatch( filename , group1 ,"Np"   , Np    , H5T_NATIVE_INT , 2 , start2 , loc_size2 , glo_size2 );
      readPatch( filename , group1 ,"Index", Index , H5T_NATIVE_INT , 2 , start2 , loc_size2 , glo_size2 );
      memcpy( theDomain->Np , Np , ((size_t) (Nr*Nz)) * sizeof(int) );

      getH5dims( filename , group2 ,"Cells", dims );
      int Nc = (int) dims[0];
      Nq = (int) dims[1];

      //Read in each radial track one at a time, because
      //you don't really know where the different radial
      //tracks belong in memory, they might not be 
      //contiguous, and they definitely are not in a rectangular block.
        for(k = 0; k < Nz; ++k)
        {
            for(j = 0; j < Nr; ++j)
            {
                int jk = j+Nr*k;
                start2[0] = Index[jk];
                start2[1] = 0;
                loc_size2[0] = Np[jk];
                loc_size2[1] = Nq;
                glo_size2[0] = Nc;
                glo_size2[1] = Nq;
                double TrackData[Np[jk]*Nq];
                readPatch(filename, group2, "Cells", TrackData,
                          H5T_NATIVE_DOUBLE, 2, start2, loc_size2, glo_size2);
                theDomain->theCells[jk] = (struct cell *) malloc(
                        ((size_t) Np[jk]) * sizeof(struct cell));
            
                for(i = 0; i < Np[jk]; ++i)
                {
                    struct cell * c = theCells[jk]+i;
                    Doub2Cell( TrackData + i*Nq , c );
                }
            }
        }

       // Setup Planets
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
                                    ((size_t) Npl) * sizeof(struct planet) );

           theDomain->pl_gas_track = (double *) malloc(
                        ((size_t) (Npl * NUM_PL_INTEGRALS)) * sizeof(double));
           theDomain->pl_kin = (double *) malloc(
                        ((size_t) (Npl * NUM_PL_KIN)) * sizeof(double));
           theDomain->pl_RK_kin = (double *) malloc(
                        ((size_t) (Npl * NUM_PL_KIN)) * sizeof(double));
           theDomain->pl_aux = (double *) malloc(
                        ((size_t) (Npl * NUM_PL_AUX)) * sizeof(double));
           theDomain->pl_RK_aux = (double *) malloc(
                        ((size_t) (Npl * NUM_PL_AUX)) * sizeof(double));
       }
      
      initializePlanets( theDomain->thePlanets );
      
      getH5dims( filename , group2 ,"Planets", dims );
      int NpDat = (int) dims[1];


      double PlanetData[Npl*NpDat];
      start2[0] = 0;
      start2[1] = 0;
      loc_size2[0] = Npl;
      loc_size2[1] = NpDat;
      glo_size2[0] = Npl;
      glo_size2[1] = NpDat;

      readPatch( filename , group2 ,"Planets", PlanetData , H5T_NATIVE_DOUBLE,
                    2, start2 , loc_size2 , glo_size2 );
      int p, q;
      for(p=0; p<Npl; ++p)
      {
          /*
           * This is a little ugly to preserve backwards-compatibility.
           * The checkpoint we're restarting from may have been made from a
           * version of the code that did not include the KIN parameters
           * or vertical position & velocity.
           */
         struct planet * pl = theDomain->thePlanets+p;
         pl->M     = PlanetData[NpDat*p + 0];
         pl->vr    = PlanetData[NpDat*p + 1];
         pl->omega = PlanetData[NpDat*p + 2];
         pl->r     = PlanetData[NpDat*p + 3];
         pl->phi   = PlanetData[NpDat*p + 4];
         pl->eps   = PlanetData[NpDat*p + 5];
         if(NpDat > 6)
            pl->type  = (int)PlanetData[NpDat*p + 6];
         else
            pl->type  = PLPOINTMASS;

         if(NpDat == 7 + NUM_PL_KIN)
         {
             pl->z = 0.0;
             pl->vz = 0.0;

             for(q=0; q<NUM_PL_KIN; q++)
                 theDomain->pl_kin[p*NUM_PL_KIN+q] = PlanetData[NpDat*p+q+7];
         }
         else if(NpDat == 9 + NUM_PL_KIN)
         {
             pl->vz = PlanetData[NpDat*p + 7];
             pl->z = PlanetData[NpDat*p + 8];
             
             for(q=0; q<NUM_PL_KIN; q++)
                 theDomain->pl_kin[p*NUM_PL_KIN+q] = PlanetData[NpDat*p+q+9];
         }
         else
         {
             pl->z = 0.0;
             pl->vz = 0.0;

             planet_init_kin(theDomain->thePlanets + p,
                             theDomain->pl_kin + p*NUM_PL_KIN);
         }
      }
      
      setPlanetsXYZ(theDomain);
      zeroAuxPlanets(theDomain);
      initializePlanetTracking(theDomain);
   }
#if USE_MPI
   MPI_Barrier(theDomain->theComm);
#endif
   }
   if( Nq != NUM_Q+NUM_FACES+1 ){ if(rank==0)printf("Ummm, I got an hdf5 read error. Check NUM_Q.\n"); exit(1); }

   
    // Setup Diagnostics & snapshot
    setup_diagnostics(theDomain);
    setup_snapshot(theDomain);

   theDomain->N_ftracks_r = get_num_rzFaces( theDomain->Nr , theDomain->Nz , 1 ); 
   theDomain->N_ftracks_z = get_num_rzFaces( theDomain->Nr , theDomain->Nz , 2 ); 
   theDomain->fIndex_r = (int *) malloc( ((size_t) (theDomain->N_ftracks_r+1))
                                        * sizeof(int) );
   theDomain->fIndex_z = (int *) malloc( ((size_t) (theDomain->N_ftracks_z+1))
                                        * sizeof(int) );

   // Reset parameters in case the restart changed things
   // This really necessary only if pointers moved around.
   setICparams( theDomain );
   setHydroParams( theDomain );
   setRiemannParams( theDomain );
   setHlldParams( theDomain );
   setOmegaParams( theDomain );
   setRotFrameParams( theDomain );
   setMetricParams( theDomain );
   setFrameParams( theDomain );
   setDiagParams( theDomain );
   setReportParams( theDomain );
   setNoiseParams( theDomain );
   setBCParams( theDomain );
   setSinkParams( theDomain );
}

