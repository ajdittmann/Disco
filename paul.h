#ifndef PAUL
#define PAUL

enum{RHO,PPP,URR,UPP,UZZ,BRR,BPP,BZZ};
enum{DDD,TAU,SRR,LLL,SZZ};

enum{PROF_TOT, PROF_DT, PROF_TIMESTEP, PROF_OUTPUT, PROF_RECON, PROF_FLUX,
     PROF_CT, PROF_SOURCE, PROF_C2P, PROF_BOUND, PROF_EXCHANGE,
     PROF_RECON_R, PROF_RECON_P, PROF_RECON_Z, 
     PROF_FLUX_R, PROF_FLUX_P, PROF_FLUX_Z,
     PROF_EXCH_NP_COUNT1, PROF_EXCH_NP_COMM1, PROF_EXCH_NP_COUNT2,
     PROF_EXCH_NP_COMM2, PROF_EXCH_NP_FIN,
     PROF_EXCH_PREP, PROF_EXCH_COMM, PROF_EXCH_FIN,
     PROF_STEP_INIT, PROF_STEP_PL, PROF_MOVE, PROF_CLEAN,
     NUM_PROF}; // NUM_PROF must be at end
enum{PLPOINTMASS, PLPW, PLSURFACEGRAV, PLSPLINE, PLWEGGC, PLQUAD,
     PLUNIFORM, PLLINEARX};
enum{COOL_NONE, COOL_BETA, COOL_BETA_RELAX};
enum{PL_M, PL_R, PL_PHI, PL_Z, PL_PR, PL_LL, PL_PZ, PL_SZ, PL_EINT};
enum{PL_SNK_M, PL_GRV_PX, PL_GRV_PY, PL_GRV_PZ, PL_GRV_JZ,
     PL_SNK_PX, PL_SNK_PY, PL_SNK_PZ, PL_SNK_JZ, PL_SNK_SZ,
     PL_SNK_MX, PL_SNK_MY, PL_SNK_MZ, PL_GRV_EGAS, PL_SNK_EGAS, PL_SNK_UGAS,
     //Up to here are all integrals
     // Everything past here are computed from the integrals above
     PL_GRV_LZ, PL_SNK_LZ, PL_GRV_K, PL_SNK_K, PL_GRV_U, PL_SNK_U, PL_SNK_EINT,
     PL_EXT_PX, PL_EXT_PY, PL_EXT_PZ, PL_EXT_JZ, PL_EXT_K, PL_EXT_U};

#if USE_MPI
#include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

// NUM_C, NUM_N, and CT_MODE are specified at compile time and defined
// with -D in the Makefile

#define NUM_Q (NUM_C+NUM_N)
#define NUM_G 2

#define NUM_PL_KIN  9
#define NUM_PL_INTEGRALS 16
#define NUM_PL_AUX (NUM_PL_INTEGRALS + 13)

//Magnetic field tracking things.  Can be set to zero if there is no MHD.
#if CT_MODE == 0        //No CT
    #define NUM_EDGES 0    
    #define NUM_FACES 0    
    #define NUM_AZ_EDGES 0 
#elif CT_MODE == 1      //2D MHD, no E^phi
    #define NUM_EDGES 4    
    #define NUM_FACES 3    
    #define NUM_AZ_EDGES 0 
#elif CT_MODE == 2      //3D MHD
    #define NUM_EDGES 8    
    #define NUM_FACES 5    
    #define NUM_AZ_EDGES 4 
#else                   //default
    #define NUM_EDGES 0    
    #define NUM_FACES 0 
    #define NUM_AZ_EDGES 0 
#endif

// Hack to touch a parameter if it is unused in a function.
#define UNUSED(x) (void)(x)

struct param_list{

   double t_min, t_max;
   int Num_R, Num_Z;
   double aspect;
   int NumRepts, NumSnaps, NumChecks;
   int Out_LogTime;

   double rmin, rmax;
   double zmin, zmax;
   double phimax;

   int NoBC_Rmin, NoBC_Rmax, NoBC_Zmin, NoBC_Zmax;

   int LogZoning, R_Periodic, Z_Periodic;
   double LogRadius;
   double MaxShort, MaxLong;
   int Mesh_Motion, Riemann_Solver, Timestep;
   int Absorb_BC, Initial_Regrid, include_atmos;

   double CFL, PLM, maxDT;
   int Cartesian_Interp;
   double Cartesian_Interp_R0;
   double Density_Floor, Pressure_Floor;

   int Exact_Mesh_Omega;
   double Exact_Mesh_Omega_Par;
   int Energy_Omega;
   double Energy_Omega_Par;
   int RotFrame;
   double RotOmega, RotD;

   double Adiabatic_Index;

   int visc_flag;
   int visc_profile;
   double viscosity;
   double visc_par;

   int isothermal_flag;
   int Cs2_Profile;
   double Cs2_Par;

   double Disk_Mach;
   double Mass_Ratio;
   double Eccentricity;
   double Drift_Rate,Drift_Exp;
   int grav2D;

   int restart_flag;
   int CT;

   int metricPar0;
   double metricPar1;
   double metricPar2;
   double metricPar3;
   double metricPar4;
   
   int initPar0;
   double initPar1;
   double initPar2;
   double initPar3;
   double initPar4;
   double initPar5;
   double initPar6;
   double initPar7;
   double initPar8;

   int noiseType;
   double noiseAbs;
   double noiseRel;

   int sinkType;
   int sinkNumber;
   double sinkPar1;
   double sinkPar2;
   double sinkPar3;
   double sinkPar4;
   double sinkPar5;
   int nozzleType;
   double nozzlePar1;
   double nozzlePar2;
   double nozzlePar3;
   double nozzlePar4;
   int coolType;
   double coolPar1;
   double coolPar2;
   double coolPar3;
   double coolPar4;

   int dampInnerType;
   int dampOuterType;
   int dampUpperType;
   int dampLowerType;
   double dampTimeInner;
   double dampLenInner;
   double dampTimeOuter;
   double dampLenOuter;
   double dampTimeLower;
   double dampLenLower;
   double dampTimeUpper;
   double dampLenUpper;

   double grav_eps;
};

struct diagnostic_avg{
   double * Qrz;
   
   double * F_r;
   double * F_z;
   double * Fvisc_r;
   double * Fvisc_z;
   double * RK_F_r;
   double * RK_F_z;
   double * RK_Fvisc_r;
   double * RK_Fvisc_z;

   double * S;
   double * Sgrav;
   double * Svisc;
   double * Ssink;
   double * Scool;
   double * Sdamp;

   double * RK_S;
   double * RK_Sgrav;
   double * RK_Svisc;
   double * RK_Ssink;
   double * RK_Scool;
   double * RK_Sdamp;

   double t_avg;
};

struct snapshot_data {
    double *Qrz;
    double *Qarr;
    int num_Qrz;
    int num_Qarr;
};

struct domain{

   struct cell ** theCells;
   struct face * theFaces_1;
   struct face * theFaces_2;
   struct planet * thePlanets;
   struct profiler *prof;
   int * Np;
   int Nr,Nz,Ng;
   int NgRa, NgRb, NgZa, NgZb;
   int N0r, N0z, Nr_glob, Nz_glob, N0r_glob, N0z_glob;
   int N_ftracks_r;
   int N_ftracks_z;
   int Npl;
   double * r_jph;
   double * z_kph;
   double phi_max;
   int * fIndex_r;
   int * fIndex_z;
   double dr0;

   time_t Wallt_init;
   int rank,size;
   int dim_rank[2];
   int dim_size[2];
   int left_rank[2];
   int right_rank[2];
#if USE_MPI
   MPI_Comm theComm;
#endif

   int planet_gas_track_synced;
   double *pl_gas_track;
   double *pl_kin;
   double *pl_RK_kin;
   double *pl_aux;
   double *pl_RK_aux;

   struct param_list theParList;
   int num_tools;
   struct diagnostic_avg theTools;
   struct snapshot_data theSnap;

   double t;
   int count_steps;
   double t_init, t_fin;
   int nrpt;
   int N_rpt;
   int nsnp;
   int N_snp;
   int nchk;
   int N_chk;

   int final_step;
   int check_plz;
   int startup;
};

struct cell{

   double prim[NUM_Q];
   double cons[NUM_Q];
   double RKcons[NUM_Q];
   double gradr[NUM_Q];
   double gradp[NUM_Q];
   double gradz[NUM_Q];
   double piph;
   double dphi;
   double wiph;
   double xyz[3];


#if NUM_EDGES > 0
   double E[NUM_EDGES];
   double B[NUM_EDGES];
#endif
#if NUM_AZ_EDGES > 0
   double E_phi[NUM_AZ_EDGES];
#endif
#if NUM_FACES > 0
   double    Phi[NUM_FACES];
   double RK_Phi[NUM_FACES];
#endif
   
   double tempDoub;

   int real;
};

struct edge{
   struct cell * LU;
   struct cell * RU;
   struct cell * LD;
   struct cell * RD;

   int Prim_Zone;
   int Alt_LR;
   int Alt_UD;

   double E_dl;
};

struct face{
   struct cell * L;
   struct cell * R;
   double dxL;
   double dxR;
   double cm[3];
   double dphi;
   double dl;
   double dA;

   double E,B;
   int LRtype;
   int flip_flag;
};

struct planet{
   double r;
   double phi; 
   double z;
   double M;
   double omega;
   double vr;
   double vz;

   double eps;

   double Uf;

   int type;

   double xyz[3];
};

struct profiler{
    clock_t ticks[NUM_PROF];
    struct timespec time[NUM_PROF];
    clock_t elapsed_ticks[NUM_PROF];
    struct timespec elapsed_time[NUM_PROF];
};

#endif
