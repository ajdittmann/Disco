#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vutil.h"
#include "rdisco.h"

#define TINY 1e-100
#define TOL 1e-10

struct Grid grid_new(char *filename)
{
    printf("New grid from %s\n", filename);
    struct Grid g = {0};

    char grid_group[256];
    char data_group[256];
    char pars_group[256];
    char opts_group[256];
    strcpy(grid_group, "Grid");
    strcpy(data_group, "Data");
    strcpy(pars_group, "Pars");
    strcpy(opts_group, "Opts");

    hsize_t dims[3];
    
    getH5dims(filename, grid_group, (char *)"r_jph", dims);
    g.Nr = dims[0] - 1;
    
    getH5dims(filename, grid_group, (char *)"z_kph", dims);
    g.Nz = dims[0] - 1;
    
    getH5dims(filename, data_group, (char *)"Cells", dims);
    g.Nq = dims[1] - 1;

    g.Np = (int *) malloc(g.Nr * g.Nz * sizeof(int));
    g.idx = (int *) malloc(g.Nr * g.Nz * sizeof(long));
    g.rjph = (double *) malloc( (g.Nr+1) * sizeof(double));
    g.zkph = (double *) malloc( (g.Nz+1) * sizeof(double));

    readSimple(filename, grid_group, (char *)"r_jph", g.rjph,
                H5T_NATIVE_DOUBLE);
    readSimple(filename, grid_group, (char *)"z_kph", g.zkph,
                H5T_NATIVE_DOUBLE);
    readSimple(filename, pars_group, (char *)"Phi_Max", &(g.pmax),
                H5T_NATIVE_DOUBLE);

    char buf[256];
    readString(filename, opts_group, (char *)"GEOMETRY", buf, 256);
    if(strlen(buf) == 0)
        strcpy(buf, "cylindrical");

    if(strcmp(buf, "cylindrical") == 0)
        g.geom = GEOM_CYL;
    else if(strcmp(buf, "cartesian") == 0)
        g.geom = GEOM_CART;
    else if(strcmp(buf, "spherical") == 0)
        g.geom = GEOM_SPH;
    else if(strcmp(buf, "prolateSpheroidal") == 0)
        g.geom = GEOM_PRLSPH;
    else
        g.geom = GEOM_CYL;


    int start[2] = {0, 0};
    int loc_size[2] = {g.Nz, g.Nr};
    int glo_size[2] = {g.Nz, g.Nr};

    readPatch(filename, grid_group, (char *)"Np", g.Np,
                H5T_NATIVE_INT, 2, start, loc_size, glo_size);
    
    long Ntot = 0;

    int jk;
    for(jk=0; jk<g.Nr*g.Nz; jk++)
        Ntot += g.Np[jk];

    g.Ntot = Ntot;

    g.piph = (double *) malloc(Ntot * sizeof(double));
    g.dphi = (double *) malloc(Ntot * sizeof(double));
    g.prim = (double *) malloc(Ntot * g.Nq * sizeof(double));

    grid_load(&g, filename);

    return g;
}

void grid_load(struct Grid *g, char filename[])
{
    printf("Loading %s\n", filename);

    char grid_group[256];
    char data_group[256];
    strcpy(grid_group, "Grid");
    strcpy(data_group, "Data");

    readSimple(filename, grid_group, (char *)"T", &(g->t), H5T_NATIVE_DOUBLE);

    int start[2] = {0, 0};
    int loc_size[2] = {g->Nz, g->Nr};
    int glo_size[2] = {g->Nz, g->Nr};

    int *glo_idx = (int *) malloc(g->Nz * g->Nr * sizeof(int));

    readPatch(filename, grid_group, (char *)"Index", glo_idx,
                H5T_NATIVE_INT, 2, start, loc_size, glo_size);

    int jk;

    for(jk=0; jk<g->Nr*g->Nz; jk++)
    {
        if(jk == 0)
            g->idx[jk] = 0;
        else
            g->idx[jk] = g->idx[jk-1] + g->Np[jk-1];

        int start_prim[2] = {glo_idx[jk], 0};
        int start_piph[2] = {glo_idx[jk], g->Nq};
        int count_prim[2] = {g->Np[jk], g->Nq};
        int count_piph[2] = {g->Np[jk], 1};
        int fdims_cells[2] = {(int)g->Ntot, g->Nq+1};

        readPatch(filename, data_group, (char *)"Cells",
                  g->prim + g->Nq * g->idx[jk], H5T_NATIVE_DOUBLE,
                  2, start_prim, count_prim, fdims_cells);
        readPatch(filename, data_group, (char *)"Cells",
                  g->piph + g->idx[jk], H5T_NATIVE_DOUBLE,
                  2, start_piph, count_piph, fdims_cells);
    }
    
    free(glo_idx);

    long i;

    for(jk=0; jk<g->Nr*g->Nz; jk++)
    {
        for(i=g->idx[jk]; i < g->idx[jk] + g->Np[jk]; i++)
        {
            long im = i > g->idx[jk] ? i - 1 : g->idx[jk] + g->Np[jk] - 1;
            g->dphi[i] = g->piph[i] - g->piph[im];
            while(g->dphi[i] > g->pmax)
                g->dphi[i] -= g->pmax;
            while(g->dphi[i] < 0.0)
                g->dphi[i] += g->pmax;
        }
    }

    for(jk=0; jk<g->Nr*g->Nz; jk++)
    {
        long i0 = g->idx[jk];
        double offset = 0;
        if(g->piph[i0] - g->dphi[i0] > 0.0)
            offset = -g->pmax;

        g->piph[i0] += offset;

        for(i=i0+1; i < i0+g->Np[jk]; i++)
        {
            if(g->piph[i] + offset < g->piph[i-1])
                offset += g->pmax;

            g->piph[i] += offset;
        }

        for(i=i0+1; i<i0+g->Np[jk]; i++)
        {
            if(g->piph[i] <= g->piph[i-1])
            {
                printf("Bad mesh fix at jk=%d:  %.3g, %.3g, ... %.3g, %.3g, ...%.3g, %.3g",
                        jk, g->piph[i0], g->piph[i0+1], g->piph[i-1],
                        g->piph[i], g->piph[g->Np[jk]-2],
                        g->piph[g->Np[jk]-1]);
                break;
            }
        }
    }
}

void grid_free(struct Grid *g)
{
    g->Nr = 0;
    g->Nz = 0;
    safe_free((void **) &(g->Np));
    safe_free((void **) &(g->idx));
    safe_free((void **) &(g->rjph));
    safe_free((void **) &(g->zkph));
    safe_free((void **) &(g->piph));
    safe_free((void **) &(g->dphi));
    safe_free((void **) &(g->prim));
}

void grid_print(struct Grid *g)
{
    printf("t = %.3g\n", g->t);
    printf("Nr = %d    Nz = %d\n", g->Nr, g->Nz);
    printf("rjph: %.3g", g->rjph[0]);
    int j;
    for(j=0; j < g->Nr; j++)
        printf(", %.3g", g->rjph[j+1]);
    printf("\n");
    printf("zkph: %.3g", g->zkph[0]);
    int k;
    for(k=0; k < g->Nz; k++)
        printf(", %.3g", g->zkph[k+1]);
    printf("\n");

    long idx = 0;

    double avg[g->Nq];
    
    int q;
    for(q=0; q<g->Nq; q++)
        avg[q] = 0.0;

    for(idx=0; idx < g->Ntot; idx++)
        for(q=0; q < g->Nq; q++)
            avg[q] += g->prim[g->Nq * idx + q];
    
    for(q=0; q<g->Nq; q++)
        avg[q] /= g->Ntot;

    printf("Prim Avg: %.3g", avg[0]);
    for(q=1; q<g->Nq; q++)
        printf(", %.3g", avg[q]);
    printf("\n");
}

void grid_pos(struct Grid *g, double xyz[3], double x[3])
{
    if(g->geom == GEOM_CYL)
        grid_pos_cyl(g, xyz, x);
}

void grid_pos_cyl(struct Grid *g, double xyz[3], double x[3])
{
    x[0] = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
    x[1] = atan2(xyz[1], xyz[0]);
    x[2] = xyz[2];
}

void grid_ray_bound(struct Grid *g, struct Ray *r, double *sin, double *sout)
{
    if(g->geom == GEOM_CYL)
        grid_ray_bound_cyl(g, r, sin, sout);
}

void grid_ray_bound_cyl(struct Grid *g, struct Ray *r,
                        double *sin, double *sout)
{
    double rmax = g->rjph[g->Nr];
    double zmin = g->zkph[0];
    double zmax = g->zkph[g->Nz];

    double nt = sqrt(r->n[0] * r->n[0] + r->n[1] * r->n[1]);
    double nz = r->n[2];

    double r_min_ray = fabs(r->x0[0] * r->n[1] - r->x0[1] * r->n[0]) / nt;
    double orig_ray = -(r->x0[0]*r->n[0] + r->x0[1]*r->n[1]) / nt;

    double r0 = sqrt(r->x0[0]*r->x0[0] + r->x0[1]*r->x0[1]);

    double sin_r, sout_r;

    if(nt < TINY || r_min_ray >= rmax)
    {
        if(r0 < rmax)
        {
            sin_r = -HUGE;
            sout_r = HUGE;
        }
        else
        {
            sin_r = HUGE;
            sout_r = -HUGE;
        }
    }
    else
    {
        sin_r = (orig_ray - sqrt((rmax-r_min_ray)*(rmax+r_min_ray))) / nt;
        sout_r = (orig_ray + sqrt((rmax-r_min_ray)*(rmax+r_min_ray))) / nt;
    }

    double sin_z, sout_z;

    if(nz > 0)
    {
        sin_z = (zmin - r->x0[2]) / nz;
        sout_z = (zmax - r->x0[2]) / nz;
    }
    else if (nz < 0)
    {
        sin_z = (zmax - r->x0[2]) / nz;
        sout_z = (zmin - r->x0[2]) / nz;
    }
    else if (r->x0[2] > zmin && r->x0[2] < zmax)
    {
        sin_z = -HUGE;
        sout_z = HUGE;
    }
    else
    {
        sin_z = HUGE;
        sout_z = -HUGE;
    }

    if(sin_r >= sout_r || sin_z >= sout_z)
    {
        *sin = HUGE;
        *sout = -HUGE;
        return;
    }

    if(sin_r < sin_z)
        *sin = sin_z;
    else
        *sin = sin_r;

    if(sout_r < sout_z)
        *sout = sout_r;
    else
        *sout = sout_z;
}

long grid_find_cell_at_pos(struct Grid *g, double xyz[3], int *k, int *j)
{
    double x[3] = {0};
    grid_pos(g, xyz, x);

    *k = search_sorted(g->zkph, g->Nz+1, x[2]) - 1;
    *j = search_sorted(g->rjph, g->Nr+1, x[0]) - 1;

    if(*k < 0 || *k >= g->Nz || *j < 0 || *j >= g->Nr)
        return -1;

    int jk = (*k) * g->Nr + (*j);

    int ia = g->idx[jk];
    int ib = g->idx[jk] + g->Np[jk]-1;

    if(x[1] < g->piph[ia] - g->dphi[ia])
        x[1] += g->pmax;
    else if(x[1] >= g->piph[ib])
        x[1] -= g->pmax;

    int i = ia + search_sorted(g->piph + ia, g->Np[jk], x[1]);

    return i;
}

struct IntersectionList intersect(struct Grid *g, struct Ray *r)
{
    struct IntersectionList l = intersectList_new(g->Nr);

    double s_min, s_max;

    grid_ray_bound(g, r, &s_min, &s_max);

    if(s_min >= s_max)
        return l;

    int i;
    double ds_avg = 0.05;
    int Nsamps = (int) ((s_max-s_min) / ds_avg + 2);
    //int Nsamps = 100;

    for(i=0; i<Nsamps; i++)
    {
        double s = s_min + (i+0.5)*(s_max-s_min) / Nsamps;
        double x[3];
        ray_xyz(r, s, x);
        int j = -1;
        int k = -1;
        long idx = grid_find_cell_at_pos(g, x, &k, &j);

        if(idx < 0 || j < 0 || k < 0
                || idx > g->Ntot || j >= g->Nr || k >= g->Nz)
            continue;

        struct Intersection hit = {.idx = idx,
                                    .k = k,
                                    .j = j,
                                    .sin = s_min+i*(s_max-s_min)/Nsamps,
                                    .sout = s_min+(i+1)*(s_max-s_min)/Nsamps};
        intersectList_push(&l, &hit);
    }

    intersectList_resize(&l, l.len);

    return l;
}

void refreshIntersections(struct Grid *g, struct Ray *rays,
                            struct IntersectionList *hit_lists, int N)
{
    printf("Refreshing Intersections\n");
    int i, j;
    for(i=0; i<N; i++)
    {
        struct Ray *r = rays + i;

        for(j=0; j<hit_lists[i].len; j++)
        {
            struct Intersection *is = hit_lists[i].buf + j;
            double s = 0.5*(is->sin + is->sout);
            double x[3];
            ray_xyz(r, s, x);
            int j = -1;
            int k = -1;
            long idx = grid_find_cell_at_pos(g, x, &k, &j);
            if(idx < 0)
                continue;
            is->idx = idx;
        }
    }
}

void render(struct Grid *g, struct Ray *rays, struct IntersectionList *ray_hits,
            float *image, int Nx, int Ny,
            void (*transferFunc)(double *, double *, int, float *))
{
    printf("Rendering\n");
    int i;
    for(i=0; i<Nx*Ny; i++)
    {
        int j;
        image[3*i] = 0;
        image[3*i+1] = 0;
        image[3*i+2] = 0;

        if(ray_hits[i].len <= 0)
            continue;

        for(j=0; j<ray_hits[i].len; j++)
        {
            long idx = ray_hits[i].buf[j].idx;
            float rgba[4];
            double s = 0.5*(ray_hits[i].buf[j].sin + ray_hits[i].buf[j].sout);
            float ds = ray_hits[i].buf[j].sout - ray_hits[i].buf[j].sin;
            //ds = 1.0f;

            double xyz[3];
            ray_xyz(rays+i, s, xyz);

            transferFunc(xyz, g->prim + g->Nq*idx, g->Nq, rgba);

            float iabs = 1.0/rgba[3];
            float tau = rgba[3] * ds;
            float srcFac = -expm1(-tau) * iabs;  // (1 - exp(-tau)) / abs
            if (tau < 1.0e-6)
                srcFac = ds;
            float transFac = exp(-tau);   // exp(-tau)

            image[3*i+0] = srcFac * rgba[0] + transFac * image[3*i+0];
            image[3*i+1] = srcFac * rgba[1] + transFac * image[3*i+1];
            image[3*i+2] = srcFac * rgba[2] + transFac * image[3*i+2];
        }
    }

}

struct IntersectionList intersectList_new(int cap)
{
    struct IntersectionList l;
    l.cap = cap;
    l.buf = (struct Intersection *) malloc(cap * sizeof(struct Intersection));
    l.len = 0;

    return l;
}

void intersectList_resize(struct IntersectionList *l, int new_cap)
{
    l->buf = (struct Intersection *)realloc(l->buf,
                new_cap * sizeof(struct Intersection));
    l->cap = new_cap;
}

void intersectList_push(struct IntersectionList *l, struct Intersection *i)
{
    while(l->len >= l->cap)
        intersectList_resize(l, 2*l->cap);

    l->buf[l->len] = *i;
    (l->len)++;
}

void intersectList_free(struct IntersectionList *l)
{
    safe_free((void **) &(l->buf));
    l->cap = 0;
    l->len = 0;
}

void ray_print(struct Ray *r)
{
    printf("x0=(%.3g, %.3g, %.3g) n=(%.3g, %.3g, %.3g)",
            r->x0[0], r->x0[1], r->x0[2], r->n[0], r->n[1], r->n[2]);
}

void ray_xyz(struct Ray *r, double s, double x[3])
{
    x[0] = r->x0[0] + s*r->n[0];
    x[1] = r->x0[1] + s*r->n[1];
    x[2] = r->x0[2] + s*r->n[2];
}

void safe_free(void **p)
{
    if(*p != NULL)
    {
        free(*p);
        *p = NULL;
    }
}

int search_sorted(double *arr, int N, double x)
{
    if(x <= arr[0])
        return 0;
    else if(x >= arr[N-1])
        return N;

    unsigned int a = 0;
    unsigned int b = N-1;
    unsigned int i = ((unsigned int) N) >> 1;

    while (b-a > 1u)
    {
        i = (b+a) >> 1;
        if (arr[i] > x)
            b = i;
        else
            a = i;
    }

    return (int) b;
}

void cross(double a[3], double b[3], double c[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

void normalize(double v[3])
{
    double ivn = 1.0/sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] *= ivn;
    v[1] *= ivn;
    v[2] *= ivn;

}

void look_at(double o[3], double target[3], double n[3])
{
    n[0] = target[0] - o[0];
    n[1] = target[1] - o[1];
    n[2] = target[2] - o[2];
    normalize(n);
}

struct Ray * ray_generate_ortho(int Nx, int Ny, double *cam_xyz, double *cam_nz,
                                double *cam_nx, double dy)
{
    struct Ray *rays = (struct Ray *) malloc(Nx*Ny * sizeof(struct Ray));

    double cam_ny[3];
    cross(cam_nz, cam_nx, cam_ny);

    double dx = Nx * dy / Ny;

    int i, j;
    for(j=0; j<Ny; j++)
    {
        double cy = (2*(j+0.5) - Ny) * dy / Ny;
        for(i=0; i<Nx; i++)
        {
            double cx = (2*(i+0.5) - Nx) * dx / Nx;

            struct Ray r = {.x0 = {cam_xyz[0] + cx*cam_nx[0] + cy*cam_ny[0],
                                    cam_xyz[1] + cx*cam_nx[1] + cy*cam_ny[1],
                                    cam_xyz[2] + cx*cam_nx[2] + cy*cam_ny[2]},
                            .n = {cam_nz[0], cam_nz[1], cam_nz[2]}};

            rays[j*Nx+i] = r;
        }
    }

    return rays;
}

struct Ray * ray_generate_persp(int Nx, int Ny, double *cam_xyz, double *cam_nz,
                                double *cam_nx, double fov_y_deg)
{
    struct Ray *rays = (struct Ray *) malloc(Nx*Ny * sizeof(struct Ray));

    double cam_ny[3];
    cross(cam_nz, cam_nx, cam_ny);

    double dy = tan(0.5 * fov_y_deg * M_PI/180.0);
    double dx = Nx * dy / Ny;

    int i, j;
    for(j=0; j<Ny; j++)
    {
        double cy = (2*(j+0.5) - Ny) * dy / Ny;
        for(i=0; i<Nx; i++)
        {
            double cx = (2*(i+0.5) - Nx) * dx / Nx;

            double n[3] = {cam_nz[0] + cx*cam_nx[0] + cy*cam_ny[0],
                            cam_nz[1] + cx*cam_nx[1] + cy*cam_ny[1],
                            cam_nz[2] + cx*cam_nx[2] + cy*cam_ny[2]};
            normalize(n);

            struct Ray r = {.x0 = {cam_xyz[0], cam_xyz[1], cam_xyz[2]},
                            .n = {n[0], n[1], n[2]}};

            rays[j*Nx+i] = r;
        }
    }

    return rays;
}
