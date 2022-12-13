
#include "paul.h"
#include "geometry.h"
#include "analysis.h"
#include "output.h"

void setup_diagnostics(struct domain *theDomain)
{
    struct diagnostic_avg *theTools = &(theDomain->theTools);
    int num_tools = num_diagnostics();

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;

    theDomain->num_tools = num_tools;
    theTools->t_avg = 0.0;

    if(num_tools == 0)
        theTools->Qrz = NULL;
    else
        theTools->Qrz = (double *) malloc(Nr * Nz * num_tools
                                          * sizeof(double));

    if(Nr == 1)
    {   
        theTools->F_r = NULL;
        theTools->Fvisc_r = NULL;
        theTools->RK_F_r = NULL;
        theTools->RK_Fvisc_r = NULL;
    }
    else
    {
        theTools->F_r = (double *) malloc((Nr-1) * Nz * NUM_Q
                                            * sizeof(double));
        theTools->Fvisc_r = (double *) malloc((Nr-1) * Nz * NUM_Q
                                            * sizeof(double));
        theTools->RK_F_r = (double *) malloc((Nr-1) * Nz * NUM_Q
                                            * sizeof(double));
        theTools->RK_Fvisc_r = (double *) malloc((Nr-1) * Nz * NUM_Q
                                            * sizeof(double));
    }
    if(Nz == 1)
    {
        theTools->F_z = NULL;
        theTools->Fvisc_z = NULL;
        theTools->RK_F_z = NULL;
        theTools->RK_Fvisc_z = NULL;
    }
    else
    {
        theTools->F_z = (double *) malloc(Nr * (Nz-1) * NUM_Q
                                                * sizeof(double));
        theTools->Fvisc_z = (double *) malloc(Nr * (Nz-1) * NUM_Q
                                                * sizeof(double));
        theTools->RK_F_z = (double *) malloc(Nr * (Nz-1) * NUM_Q
                                                * sizeof(double));
        theTools->RK_Fvisc_z = (double *) malloc(Nr * (Nz-1) * NUM_Q
                                                * sizeof(double));
    }
   
    theTools->S = (double *) malloc( Nr*Nz*NUM_Q*sizeof(double) );
    theTools->Sgrav = (double *) malloc( Nr*Nz*NUM_Q*sizeof(double) );
    theTools->Svisc = (double *) malloc( Nr*Nz*NUM_Q*sizeof(double) );
    theTools->Ssink = (double *) malloc( Nr*Nz*NUM_Q*sizeof(double) );
    theTools->Scool = (double *) malloc( Nr*Nz*NUM_Q*sizeof(double) );
    theTools->Sdamp = (double *) malloc( Nr*Nz*NUM_Q*sizeof(double) );
    theTools->RK_S = (double *) malloc( Nr*Nz*NUM_Q*sizeof(double) );
    theTools->RK_Sgrav = (double *) malloc(Nr*Nz*NUM_Q*sizeof(double));
    theTools->RK_Svisc = (double *) malloc(Nr*Nz*NUM_Q*sizeof(double));
    theTools->RK_Ssink = (double *) malloc(Nr*Nz*NUM_Q*sizeof(double));
    theTools->RK_Scool = (double *) malloc(Nr*Nz*NUM_Q*sizeof(double));
    theTools->RK_Sdamp = (double *) malloc(Nr*Nz*NUM_Q*sizeof(double));

    zero_diagnostics(theDomain);
}

void setup_snapshot(struct domain *theDomain)
{
    int num_Qrz = num_snapshot_rz();
    int num_Qarr = num_snapshot_arr();
    
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;

    struct snapshot_data *theSnap = &(theDomain->theSnap);

    theSnap->num_Qrz = num_Qrz;
    theSnap->num_Qarr = num_Qarr;

    if(num_Qrz == 0)
        theSnap->Qrz = NULL;
    else
        theSnap->Qrz = (double *) malloc(Nr * Nz * num_Qrz * sizeof(double));

    if(num_Qarr == 0)
        theSnap->Qarr = NULL;
    else
        theSnap->Qarr = (double *) malloc(Nr * Nz * num_Qarr * sizeof(double));

    zero_snapshot(theDomain);
}

void zero_diagnostics( struct domain * theDomain )
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int Nq = theDomain->num_tools;
    struct diagnostic_avg * theTools = &(theDomain->theTools);

    if(Nq > 0)
    {
        memset(theTools->Qrz, 0, Nz*Nr*Nq * sizeof(double));
    }

    if(Nr > 1)
    {
        memset(theTools->F_r,        0, Nz*(Nr-1)*NUM_Q * sizeof(double));
        memset(theTools->Fvisc_r,    0, Nz*(Nr-1)*NUM_Q * sizeof(double));
        memset(theTools->RK_F_r,     0, Nz*(Nr-1)*NUM_Q * sizeof(double));
        memset(theTools->RK_Fvisc_r, 0, Nz*(Nr-1)*NUM_Q * sizeof(double));
    }

    if(Nz > 1)
    {
        memset(theTools->F_z,        0, (Nz-1)*Nr*NUM_Q * sizeof(double));
        memset(theTools->Fvisc_z,    0, (Nz-1)*Nr*NUM_Q * sizeof(double));
        memset(theTools->RK_F_z,     0, (Nz-1)*Nr*NUM_Q * sizeof(double));
        memset(theTools->RK_Fvisc_z, 0, (Nz-1)*Nr*NUM_Q * sizeof(double));
    }


    memset(theTools->S,     0, Nz*Nr*NUM_Q * sizeof(double));
    memset(theTools->Sgrav, 0, Nz*Nr*NUM_Q * sizeof(double));
    memset(theTools->Svisc, 0, Nz*Nr*NUM_Q * sizeof(double));
    memset(theTools->Ssink, 0, Nz*Nr*NUM_Q * sizeof(double));
    memset(theTools->Scool, 0, Nz*Nr*NUM_Q * sizeof(double));
    memset(theTools->Sdamp, 0, Nz*Nr*NUM_Q * sizeof(double));

    theTools->t_avg = 0.0;
}

void zero_snapshot( struct domain * theDomain )
{
    memset(theDomain->theSnap.Qrz, 0,
            theDomain->Nr * theDomain->Nz * theDomain->theSnap.num_Qrz
            * sizeof(double));
    memset(theDomain->theSnap.Qarr, 0, theDomain->theSnap.num_Qarr
            * sizeof(double));
}

void avg_diagnostics( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int Nq = theDomain->num_tools;
   struct diagnostic_avg * theTools = &(theDomain->theTools);
   double dt = theTools->t_avg;

   int j,k,q; 
   for(k=0; k<Nz; k++){
      for( j=0 ; j<Nr ; ++j ){
         for( q=0 ; q<Nq ; ++q ){
            int iq = k*Nq*Nr + j*Nq + q; 
            theTools->Qrz[iq] /= dt; 
         }
      }
   }

   for(k=0; k<Nz; k++){
      for( j=0 ; j<Nr-1 ; ++j ){
         int jk = k*(Nr-1) + j;
         for( q=0 ; q<NUM_Q ; ++q ){
            int iq = NUM_Q*jk + q; 
            theTools->F_r[iq] /= dt; 
            theTools->Fvisc_r[iq] /= dt; 
         }
      }
   }

   for(k=0; k<Nz-1; k++){
      for( j=0 ; j<Nr ; ++j ){
         int jk = k*Nr + j;
         for( q=0 ; q<NUM_Q ; ++q ){
            int iq = NUM_Q*jk + q; 
            theTools->F_z[iq] /= dt; 
            theTools->Fvisc_z[iq] /= dt; 
         }
      }
   }

   for(k=0; k<Nz; k++){
       for(j=0; j<Nr; j++){
           int jk = k*Nr + j;
           for(q=0; q<NUM_Q; q++) {
               int iq = NUM_Q*jk + q;
               theTools->S[iq] /= dt;
               theTools->Sgrav[iq] /= dt;
               theTools->Svisc[iq] /= dt;
               theTools->Ssink[iq] /= dt;
               theTools->Scool[iq] /= dt;
               theTools->Sdamp[iq] /= dt;
           }
       }
   }
       
   theTools->t_avg = 0.0; 

}

void adjust_RK_diag(struct domain *theDomain, double RK)
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    struct diagnostic_avg * theTools = &(theDomain->theTools);
    int j,k,q; 

    for(k=0; k<Nz; k++){
        for( j=0 ; j<Nr-1 ; ++j ){
            int jk = k*(Nr-1) + j;
            for( q=0 ; q<NUM_Q ; ++q ){
                int iq = NUM_Q*jk + q; 
                theTools->F_r[iq] = (1-RK)*theTools->F_r[iq]
                                    + RK*theTools->RK_F_r[iq]; 
                theTools->Fvisc_r[iq] = (1-RK)*theTools->Fvisc_r[iq]
                                        + RK*theTools->RK_Fvisc_r[iq]; 
            }
        }
    }

    for(k=0; k<Nz-1; k++){
        for( j=0 ; j<Nr; ++j ){
            int jk = k*Nr + j;
            for( q=0 ; q<NUM_Q ; ++q ){
                int iq = NUM_Q*jk + q; 
                theTools->F_z[iq] = (1-RK)*theTools->F_z[iq]
                                    + RK*theTools->RK_F_z[iq]; 
                theTools->Fvisc_z[iq] = (1-RK)*theTools->Fvisc_z[iq]
                                        + RK*theTools->RK_Fvisc_z[iq]; 
            }
        }
    }
    for(k=0; k<Nz; k++) {
        for(j=0; j<Nr; j++) {
            int jk = k*Nr + j;
            for( q=0 ; q<NUM_Q ; ++q ){
                int iq = NUM_Q*jk + q; 
                theTools->S[iq] = (1-RK) * theTools->S[iq]
                                    + RK * theTools->RK_S[iq];
                theTools->Sgrav[iq] = (1-RK) * theTools->Sgrav[iq]
                                        + RK * theTools->RK_Sgrav[iq];
                theTools->Svisc[iq] = (1-RK) * theTools->Svisc[iq]
                                        + RK * theTools->RK_Svisc[iq];
                theTools->Ssink[iq] = (1-RK) * theTools->Ssink[iq]
                                        + RK * theTools->RK_Ssink[iq];
                theTools->Scool[iq] = (1-RK) * theTools->Scool[iq]
                                        + RK * theTools->RK_Scool[iq];
                theTools->Sdamp[iq] = (1-RK) * theTools->Sdamp[iq]
                                        + RK * theTools->RK_Sdamp[iq];
            }
        }
    }
}

void copy_RK_diag(struct domain *theDomain)
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    struct diagnostic_avg * theTools = &(theDomain->theTools);

    if(Nr > 1)
    {
        memcpy(theTools->RK_F_r, theTools->F_r,
               Nz*(Nr-1)*NUM_Q*sizeof(double));
        memcpy(theTools->RK_Fvisc_r, theTools->Fvisc_r,
               Nz*(Nr-1)*NUM_Q*sizeof(double));
    }

    if(Nz > 1)
    {
        memcpy(theTools->RK_F_z, theTools->F_z,
               (Nz-1)*Nr*NUM_Q*sizeof(double));
        memcpy(theTools->RK_Fvisc_z, theTools->Fvisc_z,
               (Nz-1)*Nr*NUM_Q*sizeof(double));
    }

    memcpy(theTools->RK_S,     theTools->S,     Nr*Nz*NUM_Q*sizeof(double));
    memcpy(theTools->RK_Sgrav, theTools->Sgrav, Nr*Nz*NUM_Q*sizeof(double));
    memcpy(theTools->RK_Svisc, theTools->Svisc, Nr*Nz*NUM_Q*sizeof(double));
    memcpy(theTools->RK_Ssink, theTools->Ssink, Nr*Nz*NUM_Q*sizeof(double));
    memcpy(theTools->RK_Scool, theTools->Scool, Nr*Nz*NUM_Q*sizeof(double));
    memcpy(theTools->RK_Sdamp, theTools->Sdamp, Nr*Nz*NUM_Q*sizeof(double));
}


void add_diagnostics( struct domain * theDomain , double dt ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Nq = theDomain->num_tools;
   struct diagnostic_avg * theTools = &(theDomain->theTools);
   struct cell ** theCells = theDomain->theCells;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   double temp_sum[Nr*Nz*Nq];
   double temp_vol[Nr*Nz];
   memset( temp_sum, 0 , Nr*Nz*Nq*sizeof(double) );
   memset( temp_vol, 0 , Nr*Nz*sizeof(double) );
   int i,j,k,q;
   int kmin = 0;
   int kmax = Nz;
   int jmin = 0;
   int jmax = Nr;

    for(k=kmin; k<kmax; k++)
    {
        for(j=jmin; j<jmax; j++)
        {
            int jk = k*Nr + j;
            for(i=0; i<Np[jk]; i++)
            {
                struct cell * c = theCells[jk]+i;
                double phip = c->piph;
                double phim = phip - c->dphi;
                double xp[3] = {r_jph[j  ] , phip , z_kph[k  ]};  
                double xm[3] = {r_jph[j-1] , phim , z_kph[k-1]};
                double xc[3];
                get_centroid_arr(xp, xm, xc);  
                double dV = get_dV(xp,xm);
                double Qrz[Nq];
                get_diagnostics( xc , c->prim , Qrz , theDomain );
                for( q=0 ; q<Nq ; ++q ) 
                    temp_sum[ jk*Nq + q ] += Qrz[q]*dV;
                temp_vol[jk] += dV;
            }
        }
    }

    for(k=kmin; k<kmax; k++)
        for(j=jmin; j<jmax; j++)
        {
            int jk = k*Nr + j;
            for(q=0; q<Nq; q++)
                theTools->Qrz[jk*Nq+q] += temp_sum[jk*Nq+q]*dt / temp_vol[jk];
        }

   theTools->t_avg += dt;
}

void snapshot(struct domain *theDomain, char filename[])
{
    if(theDomain->rank==0)
        printf("Generating Snapshot...\n");

    zero_snapshot(theDomain);
    calc_snapshot_rz(theDomain);
    calc_snapshot_arr(theDomain);

#if USE_MPI
    struct snapshot_data *theSnap = &(theDomain->theSnap);
    double *send_Q = theSnap->Qarr;

    if(theDomain->rank == 0)
        send_Q = MPI_IN_PLACE;

    MPI_Reduce(send_Q, theSnap->Qarr, theSnap->num_Qarr,
               MPI_DOUBLE, MPI_SUM, 0, theDomain->theComm);
#endif

    writeSnapshot(theDomain, filename); 
}


void calc_snapshot_rz(struct domain *theDomain)
{
    double *Qrz = theDomain->theSnap.Qrz;
    int num_Qrz = theDomain->theSnap.num_Qrz;

    if(Qrz == NULL || num_Qrz == 0)
        return;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int * Np = theDomain->Np;
   
    struct cell ** theCells = theDomain->theCells;
    double * r_jph = theDomain->r_jph;
    double * z_kph = theDomain->z_kph;

    int i,j,k,q;
    int kmin = 0;
    int kmax = Nz;
    int jmin = 0;
    int jmax = Nr;

    for(k=kmin; k<kmax; k++)
    {
        for(j=jmin; j<jmax; j++)
        {
            int jk = k*Nr + j;
                
            double Qrz_dV[num_Qrz];
            double tot_dV = 0.0;

            for(q=0; q<num_Qrz; q++)
                Qrz_dV[q] = 0.0;

            for(i=0; i<Np[jk]; i++)
            {
                struct cell * c = theCells[jk]+i;
                double phip = c->piph;
                double phim = phip - c->dphi;
                double xp[3] = {r_jph[j  ] , phip , z_kph[k  ]};  
                double xm[3] = {r_jph[j-1] , phim , z_kph[k-1]};
                
                double xc[3];
                get_centroid_arr(xp, xm, xc);  
            
                double dV = get_dV(xp,xm);
                tot_dV += dV;

                double Qrz_loc[num_Qrz];
                get_snapshot_rz(xc, c->prim, Qrz_loc, theDomain);
                for(q=0; q<num_Qrz; q++) 
                    Qrz_dV[q] += Qrz_loc[q] * dV;
            }

            double idV = 1.0 / tot_dV;
            for(q=0; q<num_Qrz; q++)
                Qrz[jk*num_Qrz + q] = Qrz_dV[q] * idV;
        }
    }
}


void calc_snapshot_arr(struct domain *theDomain)
{
    double *Qarr = theDomain->theSnap.Qarr;
    int num_Qarr = theDomain->theSnap.num_Qarr;

    if(Qarr == NULL || num_Qarr == 0)
        return;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int * Np = theDomain->Np;
   
    struct cell ** theCells = theDomain->theCells;
    double * r_jph = theDomain->r_jph;
    double * z_kph = theDomain->z_kph;

    int i,j,k,q;
    int kmin = theDomain->NgZa;
    int kmax = Nz - theDomain->NgZb;
    int jmin = theDomain->NgRa;
    int jmax = Nr - theDomain->NgRb;

    for(k=kmin; k<kmax; k++)
    {
        for(j=jmin; j<jmax; j++)
        {
            int jk = k*Nr + j;
                
            for(i=0; i<Np[jk]; i++)
            {
                struct cell * c = theCells[jk]+i;
                double phip = c->piph;
                double phim = phip - c->dphi;
                double xp[3] = {r_jph[j  ] , phip , z_kph[k  ]};  
                double xm[3] = {r_jph[j-1] , phim , z_kph[k-1]};
                
                double xc[3];
                get_centroid_arr(xp, xm, xc);  
            
                double dV = get_dV(xp,xm);

                double Qarr_loc[num_Qarr];
                get_snapshot_arr(xc, c->prim, Qarr_loc, theDomain);
                
                for(q=0; q<num_Qarr; q++) 
                    Qarr[q] += Qarr_loc[q] * dV;
            }
        }
    }
}

void free_diagnostics(struct diagnostic_avg *theTools)
{
    if(theTools->Qrz != NULL)
        free(theTools->Qrz);
    
    if(theTools->F_r != NULL)
        free(theTools->F_r);
    if(theTools->Fvisc_r != NULL)
        free(theTools->Fvisc_r);
    if(theTools->F_z != NULL)
        free(theTools->F_z);
    if(theTools->Fvisc_z != NULL)
        free(theTools->Fvisc_z);
    
    if(theTools->RK_F_r != NULL)
        free(theTools->RK_F_r);
    if(theTools->RK_Fvisc_r != NULL)
        free(theTools->RK_Fvisc_r);
    if(theTools->RK_F_z != NULL)
        free(theTools->RK_F_z);
    if(theTools->RK_Fvisc_z != NULL)
        free(theTools->RK_Fvisc_z);

    if(theTools->S != NULL)
        free(theTools->S);
    if(theTools->Sgrav != NULL)
        free(theTools->Sgrav);
    if(theTools->Svisc != NULL)
        free(theTools->Svisc);
    if(theTools->Ssink != NULL)
        free(theTools->Ssink);
    if(theTools->Scool != NULL)
        free(theTools->Scool);
    if(theTools->Sdamp != NULL)
        free(theTools->Sdamp);

    if(theTools->RK_S != NULL)
        free(theTools->RK_S);
    if(theTools->RK_Sgrav != NULL)
        free(theTools->RK_Sgrav);
    if(theTools->RK_Svisc != NULL)
        free(theTools->RK_Svisc);
    if(theTools->RK_Ssink != NULL)
        free(theTools->RK_Ssink);
    if(theTools->RK_Scool != NULL)
        free(theTools->RK_Scool);
    if(theTools->RK_Sdamp != NULL)
        free(theTools->RK_Sdamp);
}

void free_snapshot(struct snapshot_data *theSnap)
{
    if(theSnap->Qrz != NULL)
        free(theSnap->Qrz);
    if(theSnap->Qarr != NULL)
        free(theSnap->Qarr);
}


