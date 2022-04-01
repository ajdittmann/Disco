#ifndef PAUL
#define PAUL

enum{RHO,PPP,URR,UPP,UZZ,BRR,BPP,BZZ};
enum{DDD,TAU,SRR,LLL,SZZ};
enum{PLPOINTMASS, PLPW, PLSURFACEGRAV, PLSPLINE, PLWEGGC};
enum{COOL_NONE, COOL_BETA, COOL_BETA_RELAX};

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

struct param_list{

   double t_min, t_max;
   int Num_R, Num_Z;
   double aspect;
   int NumRepts, NumSnaps, NumChecks;
   int Out_LogTime;

   double rmin, rmax;
   double zmin, zmax;
   double phimax;

   int focusType;
   double focusPar1;
   double focusPar2;
   double focusPar3;
   double focusPar4;


   int NoBC_Rmin, NoBC_Rmax, NoBC_Zmin, NoBC_Zmax;

   int LogZoning, R_Periodic, Z_Periodic;
   double LogRadius;
   double MaxShort, MaxLong;
   int Mesh_Motion, Riemann_Solver, Timestep;
   int Absorb_BC, Initial_Regrid, include_atmos;

   double CFL, PLM, maxDT;
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

   int sinkShapeType;
   int sinkRateType;
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
   double t_avg;
};

struct domain{

   struct cell ** theCells;
   struct face * theFaces_1;
   struct face * theFaces_2;
   struct planet * thePlanets;
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

   time_t Wallt_init;
   int rank,size;
   int dim_rank[2];
   int dim_size[2];
   int left_rank[2];
   int right_rank[2];
#if USE_MPI
   MPI_Comm theComm;
#endif

   struct param_list theParList;
   int num_tools;
   struct diagnostic_avg theTools;

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

   double E[NUM_EDGES];
   double B[NUM_EDGES];
   double E_phi[NUM_AZ_EDGES];
   double    Phi[NUM_FACES];
   double RK_Phi[NUM_FACES];
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
   double M;
   double omega;
   double vr;
   double RK_r;
   double RK_phi;
   double RK_M;
   double RK_omega;
   double RK_vr;

   double dM;
   double RK_dM;

   double Ls;
   double RK_Ls;
   double accL;
   double RK_accL;

   double gravL;
   double RK_gravL;

   double therm;
   double RK_therm;

   double accE;
   double RK_accE;
   double gravE;
   double RK_gravE;

   double eps;
   double Fr;
   double Fp;

   int type;
};

#endif
