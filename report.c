#include "paul.h"
#include "geometry.h"
#include "planet.h"
#include "report.h"

void initializeReport(struct domain *theDomain)
{
    if( theDomain->rank != 0)
        return;

    // Open read-only to check if it exists
    FILE *rFile = fopen("report.dat", "r");    

    // if report.dat exists and we're restarting, then leave!
    if(rFile && theDomain->theParList.restart_flag )
    {
        fclose(rFile);
        return;
    }

    // Clean up.
    if(rFile)
        fclose(rFile);


    // If we're still here, then either we're restarting but report.dat
    // doesn't exist (and so we need to write the header)
    // or
    // we are not restarting, and need to write report.dat anyways.

    rFile = fopen("report.dat","w");

    //Stamp file with compiler options
    fprintf(rFile, "# GIT_VERSION %s\n", GIT_VERSION);
    fprintf(rFile, "# INITIAL %s\n", INITIAL);
    fprintf(rFile, "# HYDRO %s\n", HYDRO);
    fprintf(rFile, "# GEOMETRY %s\n", GEOMETRY);
    fprintf(rFile, "# BOUNDARY %s\n", BOUNDARY);
    fprintf(rFile, "# OUTPUT %s\n", OUTPUT);
    fprintf(rFile, "# RESTART %s\n", RESTART);
    fprintf(rFile, "# PLANETS %s\n", PLANETS);
    fprintf(rFile, "# HLLD %s\n", HLLD);
    fprintf(rFile, "# ANALYSIS %s\n", ANALYSIS);
    fprintf(rFile, "# REPORT %s\n", REPORT);
    fprintf(rFile, "# METRIC %s\n", METRIC);
    fprintf(rFile, "# FRAME %s\n", FRAME);
    fprintf(rFile, "# CT_MODE %d\n", CT_MODE);

    //Print numerical parameters helpful for parsing
    fprintf(rFile, "# NUM_C %d\n", NUM_C);
    fprintf(rFile, "# NUM_N %d\n", NUM_N);
    fprintf(rFile, "# Npl %d\n", theDomain->Npl);
    fprintf(rFile, "# NUM_PL_KIN %d\n", NUM_PL_KIN);
    fprintf(rFile, "# NUM_PL_AUX %d\n", NUM_PL_AUX);

    int N_shared = num_shared_reports();
    int N_dist_aux = num_distributed_aux_reports();
    int N_dist_int = num_distributed_integral_reports();

    fprintf(rFile, "# N_shared_reports %d\n", N_shared);
    fprintf(rFile, "# N_distributed_aux_reports %d\n", N_dist_aux);
    fprintf(rFile, "# N_distributed_integral_reports %d\n", N_dist_int);

    fclose(rFile);
}


void report( struct domain * theDomain )
{
    double t = theDomain->t;
#if USE_MPI
    MPI_Comm grid_comm = theDomain->theComm;
#endif


    double cons_tot[NUM_Q] = {0};
    double *Q_shared = NULL;
    double *Q_dist = NULL;

    //Get the number of report entries
    unsigned int N_shared = (unsigned int) num_shared_reports();
    unsigned int N_dist_aux = (unsigned int) num_distributed_aux_reports();
    unsigned int N_dist_int = (unsigned int) num_distributed_integral_reports();
    unsigned int N_dist = N_dist_aux + N_dist_int;

    //Allocate the arrays if needed
    if(N_shared > 0)
    {
        Q_shared = (double *)malloc(((size_t) N_shared) * sizeof(double));
        memset(Q_shared, 0, ((size_t) N_shared) * sizeof(double));
    }
   
    if(N_dist_aux > 0 || N_dist_int > 0)
    {
        Q_dist = (double *)malloc(((size_t) N_dist) * sizeof(double));
        memset(Q_dist, 0, ((size_t) N_dist) * sizeof(double));
    }

    // Get the shared entries
    if(N_shared > 0)
        get_shared_reports(Q_shared, theDomain);

    // Get the distributed entries that do not need to be integrated
    if(N_dist_aux > 0)
        get_distributed_aux_reports(Q_dist, theDomain);

    // Now for the integrals

    double * r_jph = theDomain->r_jph;
    double * z_kph = theDomain->z_kph;
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRa = theDomain->NgRa;
    int NgRb = theDomain->NgRb;
    int NgZa = theDomain->NgZa;
    int NgZb = theDomain->NgZb;

    struct cell ** theCells = theDomain->theCells;

    int jmin = NgRa;
    int jmax = Nr-NgRb;
    int kmin = NgZa;
    int kmax = Nz-NgZb;

    int j, k, i;

    for(k = kmin; k < kmax; k++)
    {
        for(j = jmin; j < jmax; j++)
        {
            int jk = j + Nr*k;
            for(i = 0; i < Np[jk]; i++)
            {
                struct cell * c = &(theCells[jk][i]);
            
                unsigned int q;
                for(q=0; q<NUM_Q; q++)
                    cons_tot[q] += c->cons[q];
            }
        }
    }

    if(N_dist_int > 0)
    {
        for(k = kmin; k < kmax; k++)
        {
            for(j = jmin; j < jmax; j++)
            {
                int jk = j + Nr*k;
                for(i = 0; i < Np[jk]; i++)
                {
                    struct cell * c = &(theCells[jk][i]);
            
                    double phip = c->piph;
                    double phim = phip-c->dphi;
                    double xp[3] = {r_jph[j]  ,phip,z_kph[k]  };
                    double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};

                    double x[3];
                    get_centroid_arr(xp, xm, x);
                    double dV = get_dV(xp, xm);

                    double Q[N_dist_int];

                    get_distributed_integral_reports(x, c->prim, Q, theDomain);

                    unsigned int q;
                    for(q=0; q<N_dist_int; q++)
                        Q_dist[N_dist_aux+q] += Q[q] * dV;
                }
            }
        }
    }

    int rank = theDomain->rank;

#if USE_MPI
    double *send_cons = cons_tot;
    double *send_Q = Q_dist;
    
    if(rank == 0)
    {
        send_cons = MPI_IN_PLACE;
        send_Q = MPI_IN_PLACE;
    }

    MPI_Reduce(send_cons, cons_tot, NUM_Q, MPI_DOUBLE, MPI_SUM, 0, grid_comm);
    MPI_Reduce(send_Q, Q_dist, N_dist, MPI_DOUBLE, MPI_SUM, 0, grid_comm);
#endif

    if(rank == 0)
    {
        unsigned int q;

        FILE *rFile = fopen("report.dat","a");
        fprintf(rFile, "%.15le", t);

        for(q=0; q<NUM_Q; q++)
            fprintf(rFile, " %.15le", cons_tot[q]);

        for(q=0; q<N_shared; q++)
            fprintf(rFile," %.15le", Q_shared[q]);
        
        for(q=0; q<N_dist_aux; q++)
            fprintf(rFile," %.15le", Q_dist[q]);
        
        for(q=0; q<N_dist_int; q++)
            fprintf(rFile," %.15le", Q_dist[N_dist_aux + q]);
        
        fprintf(rFile,"\n");
        fclose(rFile);
    }

    zeroAuxPlanets(theDomain);

    if(Q_shared != NULL)
        free(Q_shared);
    if(Q_dist != NULL)
        free(Q_dist);
}
