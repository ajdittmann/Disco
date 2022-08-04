enum{VAR_INT,VAR_DOUB,VAR_STR};

#include "paul.h"
#include <string.h>

int readvar( char * filename , char * varname , int vartype , void * ptr ){

   FILE * inFile = fopen( filename , "r" );
   char s[512];
   char nm[512]="";
   char s1[512];
   int found = 0;
   
   while( (fgets(s,512,inFile) != NULL) && found==0 ){
      sscanf(s,"%s ",nm);
      if( strcmp(nm,varname)==0 ){
         strcpy(s1,s);
         found=1;
      }
   }
   
   fclose( inFile );
   if( found==0 )
   {
       printf("Could not read: %s\n", varname);
       return(1);
   }

   char * s2 = s1+strlen(nm)+strspn(s1+strlen(nm),"\t :=>_");

   double temp;
   char stringval[256];

   sscanf(s2,"%lf",&temp);
   sscanf(s2,"%256s",stringval);

   if( vartype == VAR_INT ){
      *((int *)   ptr) = (int)temp;
   }else if( vartype == VAR_DOUB ){
      *((double *)ptr) = (double)temp;
   }else{
      strcpy( ptr , stringval );
   }

   return(0);
}

int read_par_file( struct domain * theDomain ){

#if USE_MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&(theDomain->rank));
   MPI_Comm_size(MPI_COMM_WORLD,&(theDomain->size));
#else
   theDomain->rank = 0;
   theDomain->size = 1;
#endif

   int rank = theDomain->rank;
   int size = theDomain->size;

   struct param_list * theList = &( theDomain->theParList );

   char pfile[] = "in.par";

   int err=0;  

   int tTimes_2pi;
   int pTimes_2pi;
   int zTimes_pi = 0;

   int nrank;
   for( nrank=0 ; nrank<size ; ++nrank ){
      if( rank==nrank ){
         err += readvar( pfile , "Num_R"                 , VAR_INT  , &(theList->Num_R)           );
         err += readvar( pfile , "Num_Z"                 , VAR_INT  , &(theList->Num_Z)           );
         err += readvar( pfile , "aspect"                , VAR_DOUB , &(theList->aspect)          );
         err += readvar( pfile , "Num_Reports"           , VAR_INT  , &(theList->NumRepts)        );
         err += readvar( pfile , "Num_Snapshots"         , VAR_INT  , &(theList->NumSnaps)        );
         err += readvar( pfile , "Num_Checkpoints"       , VAR_INT  , &(theList->NumChecks)       );
         err += readvar( pfile , "T_Start"               , VAR_DOUB , &(theList->t_min)           );
         err += readvar( pfile , "T_End"                 , VAR_DOUB , &(theList->t_max)           );
         err += readvar( pfile , "R_Min"                 , VAR_DOUB , &(theList->rmin)            );
         err += readvar( pfile , "R_Max"                 , VAR_DOUB , &(theList->rmax)            );
         err += readvar( pfile , "Z_Min"                 , VAR_DOUB , &(theList->zmin)            );
         err += readvar( pfile , "Z_Max"                 , VAR_DOUB , &(theList->zmax)            );
         err += readvar( pfile , "Phi_Max"               , VAR_DOUB , &(theList->phimax)          );
         err += readvar( pfile , "Use_Logtime"           , VAR_INT  , &(theList->Out_LogTime)     );
         err += readvar( pfile , "Log_Zoning"            , VAR_INT  , &(theList->LogZoning)       );
         err += readvar( pfile , "R_Periodic"            , VAR_INT  , &(theList->R_Periodic)      );
         err += readvar( pfile , "Z_Periodic"            , VAR_INT  , &(theList->Z_Periodic)      );
         err += readvar(pfile, "Focus_Type", VAR_INT , &(theList->focusType));
         err += readvar(pfile, "Focus_Par1", VAR_DOUB, &(theList->focusPar1));
         err += readvar(pfile, "Focus_Par2", VAR_DOUB, &(theList->focusPar2));
         err += readvar(pfile, "Focus_Par3", VAR_DOUB, &(theList->focusPar3));
         err += readvar(pfile, "Focus_Par4", VAR_DOUB, &(theList->focusPar4));
         err += readvar( pfile , "NoBC_Rmin"             , VAR_INT , &(theList->NoBC_Rmin)       );
         err += readvar( pfile , "NoBC_Rmax"             , VAR_INT , &(theList->NoBC_Rmax)       );
         err += readvar( pfile , "NoBC_Zmin"             , VAR_INT , &(theList->NoBC_Zmin)       );
         err += readvar( pfile , "NoBC_Zmax"             , VAR_INT , &(theList->NoBC_Zmax)       );
         err += readvar( pfile , "Log_Radius"            , VAR_DOUB , &(theList->LogRadius)       );
         err += readvar( pfile , "Max_Aspect_Short"      , VAR_DOUB , &(theList->MaxShort)        );
         err += readvar( pfile , "Max_Aspect_Long"       , VAR_DOUB , &(theList->MaxLong)         );
         err += readvar( pfile , "Timestep"                   , VAR_INT , &(theList->Timestep)             );
         err += readvar( pfile , "CFL"                   , VAR_DOUB , &(theList->CFL)             );
         err += readvar( pfile , "Max_DT"                   , VAR_DOUB , &(theList->maxDT)             );
         err += readvar( pfile , "PLM"                   , VAR_DOUB , &(theList->PLM)             );
         err += readvar( pfile , "Cartesian_Interp"     , VAR_INT , &(theList->Cartesian_Interp)             );
         err += readvar( pfile , "Cartesian_Interp_R0"  , VAR_DOUB , &(theList->Cartesian_Interp_R0)             );
         err += readvar( pfile , "Adiabatic_Index"       , VAR_DOUB , &(theList->Adiabatic_Index) );
         err += readvar( pfile , "Isothermal"            , VAR_INT  , &(theList->isothermal_flag) );
         err += readvar( pfile , "Cs2_Profile"      , VAR_INT  , &(theList->Cs2_Profile)     );
         err += readvar( pfile , "Cs2_Par"  , VAR_DOUB  , &(theList->Cs2_Par)     );
         err += readvar( pfile , "Density_Floor"         , VAR_DOUB , &(theList->Density_Floor)   );
         err += readvar( pfile , "Pressure_Floor"        , VAR_DOUB , &(theList->Pressure_Floor)  );
         err += readvar( pfile , "Mesh_Motion"           , VAR_INT  , &(theList->Mesh_Motion)     );
         err += readvar( pfile , "Exact_Mesh_Omega"      , VAR_INT  , &(theList->Exact_Mesh_Omega)     );
         err += readvar( pfile , "Exact_Mesh_Omega_Par"  , VAR_DOUB  , &(theList->Exact_Mesh_Omega_Par)     );
         err += readvar( pfile , "Energy_Omega"      , VAR_INT  , &(theList->Energy_Omega)     );
         err += readvar( pfile , "Energy_Omega_Par"  , VAR_DOUB  , &(theList->Energy_Omega_Par)     );
         err += readvar( pfile , "RotFrame"             , VAR_INT  , &(theList->RotFrame)     );
         err += readvar( pfile , "RotOmega"             , VAR_DOUB  , &(theList->RotOmega)     );
         err += readvar( pfile , "RotD"                 , VAR_DOUB  , &(theList->RotD)     );
         err += readvar( pfile , "Riemann_Solver"        , VAR_INT  , &(theList->Riemann_Solver)  );
         err += readvar( pfile , "Initial_Regrid"        , VAR_INT  , &(theList->Initial_Regrid)  );
         err += readvar( pfile , "Restart"               , VAR_INT  , &(theList->restart_flag)    );
         err += readvar( pfile , "Use_Viscosity"         , VAR_INT  , &(theList->visc_flag)       );
         err += readvar( pfile , "Viscosity"             , VAR_DOUB , &(theList->viscosity)       );
         err += readvar( pfile , "Viscosity_Profile"          , VAR_INT  , &(theList->visc_profile)      );
         err += readvar( pfile , "Viscosity_Par"          , VAR_DOUB  , &(theList->visc_par)      );
         err += readvar( pfile , "Include_Atmos"         , VAR_INT  , &(theList->include_atmos)   );
         err += readvar( pfile , "T_Times_2pi"           , VAR_INT  , &tTimes_2pi );
         err += readvar( pfile , "P_Times_2pi"           , VAR_INT  , &pTimes_2pi );
         err += readvar( pfile , "Z_Times_pi"           , VAR_INT  , &zTimes_pi );
         err += readvar( pfile , "Mach_Number"           , VAR_DOUB , &(theList->Disk_Mach)       );
         err += readvar( pfile , "Mass_Ratio"            , VAR_DOUB , &(theList->Mass_Ratio)      );
         err += readvar( pfile , "Eccentricity"          , VAR_DOUB , &(theList->Eccentricity)    );
         err += readvar( pfile , "Drift_Rate"            , VAR_DOUB , &(theList->Drift_Rate)      );
         err += readvar( pfile , "Drift_Exp"             , VAR_DOUB , &(theList->Drift_Exp)       );
         err += readvar( pfile , "Grav2D"                , VAR_INT  , &(theList->grav2D)       );
         err += readvar( pfile , "Constrained_Transport" , VAR_INT  , &(theList->CT)              );
         err += readvar( pfile , "Metric_Par0" , VAR_INT  , &(theList->metricPar0));
         err += readvar( pfile , "Metric_Par1" , VAR_DOUB  , &(theList->metricPar1));
         err += readvar( pfile , "Metric_Par2" , VAR_DOUB  , &(theList->metricPar2));
         err += readvar( pfile , "Metric_Par3" , VAR_DOUB  , &(theList->metricPar3));
         err += readvar( pfile , "Metric_Par4" , VAR_DOUB  , &(theList->metricPar4));
         err += readvar( pfile , "Init_Par0" , VAR_INT  , &(theList->initPar0));
         err += readvar( pfile , "Init_Par1" , VAR_DOUB  , &(theList->initPar1));
         err += readvar( pfile , "Init_Par2" , VAR_DOUB  , &(theList->initPar2));
         err += readvar( pfile , "Init_Par3" , VAR_DOUB  , &(theList->initPar3));
         err += readvar( pfile , "Init_Par4" , VAR_DOUB  , &(theList->initPar4));
         err += readvar( pfile , "Init_Par5" , VAR_DOUB  , &(theList->initPar5));
         err += readvar( pfile , "Init_Par6" , VAR_DOUB  , &(theList->initPar6));
         err += readvar( pfile , "Init_Par7" , VAR_DOUB  , &(theList->initPar7));
         err += readvar( pfile , "Init_Par8" , VAR_DOUB  , &(theList->initPar8));
         err += readvar( pfile , "Noise_Type" , VAR_INT  , &(theList->noiseType));
         err += readvar( pfile , "Noise_Abs" , VAR_DOUB  , &(theList->noiseAbs));
         err += readvar( pfile , "Noise_Rel" , VAR_DOUB  , &(theList->noiseRel));
         err += readvar(pfile, "Sink_Shape_Type", VAR_INT , &(theList->sinkShapeType));
         err += readvar(pfile, "Sink_Rate_Type", VAR_INT , &(theList->sinkRateType));
         err += readvar(pfile, "Sink_Number", VAR_INT , &(theList->sinkNumber));
         err += readvar(pfile, "Sink_Par1", VAR_DOUB, &(theList->sinkPar1));
         err += readvar(pfile, "Sink_Par2", VAR_DOUB, &(theList->sinkPar2));
         err += readvar(pfile, "Sink_Par3", VAR_DOUB, &(theList->sinkPar3));
         err += readvar(pfile, "Sink_Par4", VAR_DOUB, &(theList->sinkPar4));
         err += readvar(pfile, "Sink_Par5", VAR_DOUB, &(theList->sinkPar5));
         err += readvar(pfile, "Nozzle_Type", VAR_INT , &(theList->nozzleType));
         err += readvar(pfile, "Nozzle_Par1", VAR_DOUB, &(theList->nozzlePar1));
         err += readvar(pfile, "Nozzle_Par2", VAR_DOUB, &(theList->nozzlePar2));
         err += readvar(pfile, "Nozzle_Par3", VAR_DOUB, &(theList->nozzlePar3));
         err += readvar(pfile, "Nozzle_Par4", VAR_DOUB, &(theList->nozzlePar4));
         err += readvar(pfile, "Cool_Type", VAR_INT , &(theList->coolType));
         err += readvar(pfile, "Cool_Par1", VAR_DOUB, &(theList->coolPar1));
         err += readvar(pfile, "Cool_Par2", VAR_DOUB, &(theList->coolPar2));
         err += readvar(pfile, "Cool_Par3", VAR_DOUB, &(theList->coolPar3));
         err += readvar(pfile, "Cool_Par4", VAR_DOUB, &(theList->coolPar4));
         err += readvar(pfile, "DampInner_Type", VAR_INT , &(theList->dampInnerType));
         err += readvar(pfile, "DampInner_Time", VAR_DOUB, &(theList->dampTimeInner));
         err += readvar(pfile, "DampInner_Len", VAR_DOUB, &(theList->dampLenInner));
         err += readvar(pfile, "DampOuter_Type", VAR_INT , &(theList->dampOuterType));
         err += readvar(pfile, "DampOuter_Time", VAR_DOUB, &(theList->dampTimeOuter));
         err += readvar(pfile, "DampOuter_Len", VAR_DOUB, &(theList->dampLenOuter));
         err += readvar(pfile, "DampUpper_Type", VAR_INT , &(theList->dampUpperType));
         err += readvar(pfile, "DampUpper_Time", VAR_DOUB, &(theList->dampTimeUpper));
         err += readvar(pfile, "DampUpper_Len", VAR_DOUB, &(theList->dampLenUpper));
         err += readvar(pfile, "DampLower_Type", VAR_INT , &(theList->dampLowerType));
         err += readvar(pfile, "DampLower_Time", VAR_DOUB, &(theList->dampTimeLower));
         err += readvar(pfile, "DampLower_Len", VAR_DOUB, &(theList->dampLenLower));
         err += readvar(pfile, "Softening", VAR_DOUB, &(theList->grav_eps));



      }
#if USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
   }

   if( tTimes_2pi ){
      theList->t_min *= 2.*M_PI;
      theList->t_max *= 2.*M_PI;
   }
   if( pTimes_2pi ){
      theList->phimax *= 2.*M_PI;
   }
   if( zTimes_pi ){
      theList->zmin *= M_PI;
      theList->zmax *= M_PI;
   }

   int errtot;
#if USE_MPI
   MPI_Allreduce( &err , &errtot , 1 , MPI_INT , MPI_SUM , MPI_COMM_WORLD );
#else
   errtot = err;
#endif

   if( errtot > 0 ){
      printf("Read Failed\n");
      return(1);
   }

   return(0);

}


