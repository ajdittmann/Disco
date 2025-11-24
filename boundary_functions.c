
#include "disco.h"
#include <string.h>
#include "geometry.h"
#include "hydro.h"
#include "boundary.h"

#define R_HOR 1.0

static double phi_max = 0.0;

void initial( double * , double * );
void subtract_omega( double * );

void setBCParams(struct domain *theDomain)
{
    phi_max = theDomain->phi_max;
}


void set_cell_init(struct cell *c, const double *r_jph, const double *z_kph,
                   int j, int k)
{
    double xm[3] = {r_jph[j-1], c->piph - c->dphi, z_kph[k-1]};
    double xp[3] = {r_jph[j  ], c->piph          , z_kph[k  ]};
    double x[3];
    get_centroid_arr(xp, xm, x);

    initial(c->prim, x);
    subtract_omega(c->prim);
}

void set_cell_init_q(struct cell *c, const double *r_jph, const double *z_kph, 
                        int j, int k, int *qarr, int nq)
{
    double r = get_centroid(r_jph[j], r_jph[j-1], 1);
    double z = get_centroid(z_kph[k], z_kph[k-1], 2);
    double phi = c->piph - 0.5*c->dphi;
    double x[3] = {r, phi, z};
    double temp_prim[NUM_Q];
    initial(temp_prim, x);
    subtract_omega(temp_prim);

    int iq;
    for(iq=0; iq<nq; iq++)
        c->prim[qarr[iq]] = temp_prim[qarr[iq]];
}

void set_cells_copy_distant(struct cell *c, int Np, struct cell *c1, int Np1)
{
    // Copies annulus c1 to c.

    int i, q;

    // Clear annulus.
    for(i=0; i<Np; i++)
        for(q=0; q<NUM_Q; q++)
            c[i].prim[q] = 0.0;

    // Find first cell in c1 that intersects with c[0].
    int i1 = 0;
    double phi0 = c[0].piph - c[0].dphi;
    for(i1=0; i1 < Np1; i1++)
    {
        double dphip = get_signed_dp(c1[i1].piph, phi0);
        double dphim = dphip - c1[i1].dphi;

        if(dphip > 0 && dphim <= 0.0)
            break;
    }

    if(i1 >= Np1)
        i1 = 0;

    //Loop through all cells c, add contributions from all neighbours in c1.
    for(i=0; i<Np; i++)
    {
        double phip = c[i].piph;
        double phim = c[i].piph - c[i].dphi;
        double phip1, phim1;
        do
        {
            phip1 = c1[i1].piph;
            phim1 = c1[i1].piph - c1[i1].dphi;

            while(phip1 < phim)
            {
                phip1 += phi_max;
                phim1 += phi_max;
            }
            while(phim1 > phip)
            {
                phip1 -= phi_max;
                phim1 -= phi_max;
            }
            double phi1 = phip < phip1 ? phip : phip1;
            double phi2 = phim > phim1 ? phim : phim1;
            double dphi = phi1-phi2;

            if(dphi< 0.0)
                printf("WHOA %d %d %.12lg\n", i, i1, dphi);
            for(q=0; q<NUM_Q; q++)
                c[i].prim[q] += c1[i1].prim[q] * dphi;

            i1 = i1 < Np1-1 ? i1+1 : 0;
            if(i1 == Np1)
                i1 = 0;
        }
        while(phip1 < phip);
        i1 = i1 > 0 ? i1-1 : Np1-1;

    }

    // Divide by total area.
    for(i=0; i<Np; i++)
        for(q=0; q<NUM_Q; q++)
            c[i].prim[q] /= c[i].dphi;
}

void boundary_fixed_rinn( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRa = theDomain->NgRa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[0] == 0 )
    {
        for(k=0; k<Nz; k++)
            for(j=0; j<NgRa; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
            }  
    }
}

void boundary_fixed_rout( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRb = theDomain->NgRb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(k=0; k<Nz; k++)
            for(j=Nr-NgRb; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
            }  
    }
}

void boundary_fixed_zbot( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;
    int NgZa = theDomain->NgZa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[1] == 0)
    {
        for(k=0; k<NgZa; k++)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
            }  
    }
}
void boundary_fixed_ztop( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgZb = theDomain->NgZb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(k=Nz-NgZb; k<Nz; k++)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
            }  
    }
}

void boundary_zerograd_rinn( struct domain *theDomain, int diode)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRa = theDomain->NgRa;

    int *dim_rank = theDomain->dim_rank;

    int j,k;

    if(dim_rank[0] == 0 )
    {
        for(k=0; k<Nz; k++)
            for(j=NgRa-1; j>=0; j--)
            {
                int jk = j+Nr*k;
                int jk1 = j+1+Nr*k;
                set_cells_copy_distant(theCells[jk], Np[jk],
                                        theCells[jk1], Np[jk1]);

                if(diode)
                {
                    int i;
                    for(i=0; i<Np[jk]; i++)
                    {
                        if(theCells[jk][i].prim[URR] > 0)
                            theCells[jk][i].prim[URR] = 0.0;
                    }
                }
            }
    }
}

void boundary_zerograd_rout( struct domain *theDomain, int diode)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRb = theDomain->NgRb;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(k=0; k<Nz; k++)
            for(j=Nr-NgRb; j<Nr; j++)
            {
                int jk = j + Nr*k;
                int jk1 = j-1+Nr*k;
                set_cells_copy_distant(theCells[jk], Np[jk],
                                        theCells[jk1], Np[jk1]);

                if(diode)
                {
                    int i;
                    for(i=0; i<Np[jk]; i++)
                    {
                        if(theCells[jk][i].prim[URR] < 0)
                            theCells[jk][i].prim[URR] = 0.0;
                    }
                }
            }
    }
}

void boundary_zerograd_zbot( struct domain *theDomain, int diode)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;
    int NgZa = theDomain->NgZa;

    int *dim_rank = theDomain->dim_rank;

    int j,k;

    if(dim_rank[1] == 0)
    {
        for(j=0; j<Nr; j++)
            for(k=NgZa-1; k>=0; k--)
            {
                int jk = j+Nr*k;
                int jk1 = j+Nr*(k+1);
                set_cells_copy_distant(theCells[jk], Np[jk],
                                        theCells[jk1], Np[jk1]);
                
                if(diode)
                {
                    int i;
                    for(i=0; i<Np[jk]; i++)
                    {
                        if(theCells[jk][i].prim[UZZ] > 0)
                            theCells[jk][i].prim[UZZ] = 0.0;
                    }
                }
            }
    }
}

void boundary_zerograd_ztop( struct domain *theDomain, int diode)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgZb = theDomain->NgZb;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(j=0; j<Nr; j++)
            for(k=Nz-NgZb; k<Nz; k++)
            {
                int jk = j+Nr*k;
                int jk1 = j+Nr*(k-1);
                set_cells_copy_distant(theCells[jk], Np[jk],
                                        theCells[jk1], Np[jk1]);
                
                if(diode)
                {
                    int i;
                    for(i=0; i<Np[jk]; i++)
                    {
                        if(theCells[jk][i].prim[UZZ] < 0)
                            theCells[jk][i].prim[UZZ] = 0.0;
                    }
                }
            }
    }
}

void boundary_reflect_rinn( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRa = theDomain->NgRa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[0] == 0 )
    {
        for(k=0; k<Nz; k++)
        {
            double z = get_centroid(z_kph[k], z_kph[k-1], 2);
            for(j=NgRa-1; j>=0; j--)
            {
                int jk = j+Nr*k;
                
                int j1 = 2*NgRa-j-1;
                int jk1 = j1 + Nr*k;
                
                set_cells_copy_distant(theCells[jk], Np[jk], 
                                        theCells[jk1], Np[jk1]);
                    
                double r = get_centroid(r_jph[j], r_jph[j-1], 1);

                for(i=0; i<Np[jk]; i++)
                {
                    struct cell *c = &(theCells[jk][i]);
                    double phi = c->piph - 0.5*c->dphi;
                    double x[3] = {r, phi, z};
                    reflect_prims(theCells[jk][i].prim, x, 0);
                }
            }
        }
    }
}

void boundary_reflect_rout( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRb = theDomain->NgRb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(k=0; k<Nz; k++)
        {
            double z = get_centroid(z_kph[k], z_kph[k-1], 2);

            for(j=Nr-NgRb; j<Nr; j++)
            {
                int jk = j+Nr*k;
                
                int j1 = 2*(Nr-NgRb) - j - 1;
                int jk1 = j1 + Nr*k;
                
                set_cells_copy_distant(theCells[jk], Np[jk], 
                                        theCells[jk1], Np[jk1]);
                
                double r = get_centroid(r_jph[j], r_jph[j-1], 1);

                for(i=0; i<Np[jk]; i++)
                {
                    struct cell *c = &(theCells[jk][i]);
                    double phi = c->piph - 0.5*c->dphi;
                    double x[3] = {r, phi, z};
                    reflect_prims(theCells[jk][i].prim, x, 0);
                }
            }
        }
    }
}

void boundary_reflect_zbot( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;
    int NgZa = theDomain->NgZa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[1] == 0)
    {
        for(k=NgZa-1; k>=0; k--)
        {
            double z = get_centroid(z_kph[k], z_kph[k-1], 2);

            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                
                int k1 = 2*NgZa-k-1;
                int jk1 = j + Nr*k1;
                
                set_cells_copy_distant(theCells[jk], Np[jk], 
                                        theCells[jk1], Np[jk1]);

                double r = get_centroid(r_jph[j], r_jph[j-1], 1);

                for(i=0; i<Np[jk]; i++)
                {
                    struct cell *c = &(theCells[jk][i]);
                    double phi = c->piph - 0.5*c->dphi;
                    double x[3] = {r, phi, z};
                    reflect_prims(theCells[jk][i].prim, x, 2);
                }
            }
        }
    }
}

void boundary_reflect_ztop( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgZb = theDomain->NgZb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(k=Nz-NgZb; k<Nz; k++)
        {
            double z = get_centroid(z_kph[k], z_kph[k-1], 2);

            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                
                int k1 = 2*(Nz-NgZb) - k - 1;
                int jk1 = j + Nr*k1;
                
                set_cells_copy_distant(theCells[jk], Np[jk], 
                                        theCells[jk1], Np[jk1]);

                double r = get_centroid(r_jph[j], r_jph[j-1], 1);

                for(i=0; i<Np[jk]; i++)
                {
                    struct cell *c = &(theCells[jk][i]);
                    double phi = c->piph - 0.5*c->dphi;
                    double x[3] = {r, phi, z};
                    reflect_prims(theCells[jk][i].prim, x, 2);
                }
            }
        }
    }
}

void boundary_fixed_horizon( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int i,j,k;

    for(k=0; k<Nz; k++)
    {
        double zm = z_kph[k-1];
        double zp = z_kph[k];
        double zo = fabs(zp)>fabs(zm) ? zp : zm;

        if(fabs(zo) > R_HOR)
            continue;

        for(j=0; j<Nr; j++)
        {
            double ro = r_jph[j];

            double R = sqrt(zo*zo + ro*ro);
            
            if(R < R_HOR)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                {
                    set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
                    theCells[jk][i].real = 0;
                }
            }  
        }
    }
}

void boundary_fixed_q_rinn( struct domain *theDomain, int *q, int nq)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRa = theDomain->NgRa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[0] == 0 )
    {
        for(k=0; k<Nz; k++)
            for(j=0; j<NgRa; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init_q(&(theCells[jk][i]), r_jph, z_kph, j, k,
                                    q, nq);
            }  
    }
}

void boundary_fixed_q_rout( struct domain *theDomain, int *q, int nq)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRb = theDomain->NgRb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(k=0; k<Nz; k++)
            for(j=Nr-NgRb; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init_q(&(theCells[jk][i]), r_jph, z_kph, j, k,
                                    q, nq);
            }  
    }
}

void boundary_fixed_q_zbot( struct domain *theDomain, int *q, int nq)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;
    int NgZa = theDomain->NgZa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[1] == 0)
    {
        for(k=0; k<NgZa; k++)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init_q(&(theCells[jk][i]), r_jph, z_kph, j, k,
                                    q, nq);
            }  
    }
}
void boundary_fixed_q_ztop( struct domain *theDomain, int *q, int nq)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgZb = theDomain->NgZb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(k=Nz-NgZb; k<Nz; k++)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init_q(&(theCells[jk][i]), r_jph, z_kph, j, k, 
                                    q, nq);
            }  
    }
}

void boundary_fixed_phi_rinn( struct domain *theDomain, double phia,
                              double phib)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRa = theDomain->NgRa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[0] == 0 )
    {
        for(k=0; k<Nz; k++)
            for(j=0; j<NgRa; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                {
                    double phi = theCells[jk][i].piph
                                    - 0.5*theCells[jk][i].dphi;
                    while(phi < phia)
                        phi += phi_max;
                    if(phi < phib)
                        set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
                }
            }  
    }
}

void boundary_fixed_phi_rout( struct domain *theDomain, double phia,
                            double phib)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRb = theDomain->NgRb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(k=0; k<Nz; k++)
            for(j=Nr-NgRb; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                {
                    double phi = theCells[jk][i].piph
                                    - 0.5*theCells[jk][i].dphi;
                    while(phi < phia)
                        phi += phi_max;
                    if(phi < phib)
                        set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
                }
            }  
    }
}

void boundary_fixed_phi_zbot( struct domain *theDomain, double phia,
                              double phib)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;
    int NgZa = theDomain->NgZa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[1] == 0)
    {
        for(k=0; k<NgZa; k++)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                {
                    double phi = theCells[jk][i].piph
                                    - 0.5*theCells[jk][i].dphi;
                    while(phi < phia)
                        phi += phi_max;
                    if(phi < phib)
                        set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
                }
            }  
    }
}
void boundary_fixed_phi_ztop( struct domain *theDomain, double phia,
                              double phib)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgZb = theDomain->NgZb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(k=Nz-NgZb; k<Nz; k++)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                {
                    double phi = theCells[jk][i].piph
                                    - 0.5*theCells[jk][i].dphi;
                    while(phi < phia)
                        phi += phi_max;
                    if(phi < phib)
                        set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
                }
            }  
    }
}

void boundary_noslip_rinn( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRa = theDomain->NgRa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[0] == 0 )
    {
        for(k=0; k<Nz; k++)
        {
            double z = get_centroid(z_kph[k], z_kph[k-1], 2);
            for(j=NgRa-1; j>=0; j--)
            {
                int jk = j+Nr*k;
                
                int j1 = 2*NgRa-j-1;
                int jk1 = j1 + Nr*k;
                
                set_cells_copy_distant(theCells[jk], Np[jk], 
                                        theCells[jk1], Np[jk1]);
                    
                double r = get_centroid(r_jph[j], r_jph[j-1], 1);

                for(i=0; i<Np[jk]; i++)
                {
                    struct cell *c = &(theCells[jk][i]);
                    double phi = c->piph - 0.5*c->dphi;
                    double x[3] = {r, phi, z};
                    reflect_prims(theCells[jk][i].prim, x, 0);
                    reflect_prims(theCells[jk][i].prim, x, 1);
                    reflect_prims(theCells[jk][i].prim, x, 2);
                }
            }
        }
    }
}

void boundary_noslip_rout( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRb = theDomain->NgRb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(k=0; k<Nz; k++)
        {
            double z = get_centroid(z_kph[k], z_kph[k-1], 2);

            for(j=Nr-NgRb; j<Nr; j++)
            {
                int jk = j+Nr*k;
                
                int j1 = 2*(Nr-NgRb) - j - 1;
                int jk1 = j1 + Nr*k;
                
                set_cells_copy_distant(theCells[jk], Np[jk], 
                                        theCells[jk1], Np[jk1]);
                
                double r = get_centroid(r_jph[j], r_jph[j-1], 1);

                for(i=0; i<Np[jk]; i++)
                {
                    struct cell *c = &(theCells[jk][i]);
                    double phi = c->piph - 0.5*c->dphi;
                    double x[3] = {r, phi, z};
                    reflect_prims(theCells[jk][i].prim, x, 0);
                    reflect_prims(theCells[jk][i].prim, x, 1);
                    reflect_prims(theCells[jk][i].prim, x, 2);
                }
            }
        }
    }
}

void boundary_noslip_zbot( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;
    int NgZa = theDomain->NgZa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[1] == 0)
    {
        for(k=NgZa-1; k>=0; k--)
        {
            double z = get_centroid(z_kph[k], z_kph[k-1], 2);

            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                
                int k1 = 2*NgZa-k-1;
                int jk1 = j + Nr*k1;
                
                set_cells_copy_distant(theCells[jk], Np[jk], 
                                        theCells[jk1], Np[jk1]);

                double r = get_centroid(r_jph[j], r_jph[j-1], 1);

                for(i=0; i<Np[jk]; i++)
                {
                    struct cell *c = &(theCells[jk][i]);
                    double phi = c->piph - 0.5*c->dphi;
                    double x[3] = {r, phi, z};
                    reflect_prims(theCells[jk][i].prim, x, 0);
                    reflect_prims(theCells[jk][i].prim, x, 1);
                    reflect_prims(theCells[jk][i].prim, x, 2);
                }
            }
        }
    }
}

void boundary_noslip_ztop( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgZb = theDomain->NgZb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(k=Nz-NgZb; k<Nz; k++)
        {
            double z = get_centroid(z_kph[k], z_kph[k-1], 2);

            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                
                int k1 = 2*(Nz-NgZb) - k - 1;
                int jk1 = j + Nr*k1;
                
                set_cells_copy_distant(theCells[jk], Np[jk], 
                                        theCells[jk1], Np[jk1]);

                double r = get_centroid(r_jph[j], r_jph[j-1], 1);

                for(i=0; i<Np[jk]; i++)
                {
                    struct cell *c = &(theCells[jk][i]);
                    double phi = c->piph - 0.5*c->dphi;
                    double x[3] = {r, phi, z};
                    reflect_prims(theCells[jk][i].prim, x, 0);
                    reflect_prims(theCells[jk][i].prim, x, 1);
                    reflect_prims(theCells[jk][i].prim, x, 2);
                }
            }
        }
    }
}

