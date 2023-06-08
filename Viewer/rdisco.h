#ifndef RDISCO_H
#define RDISCO_H

enum{GEOM_CART, GEOM_CYL, GEOM_SPH, GEOM_PRLSPH};

struct Ray {
    double x0[3];
    double n[3];
};

struct Intersection {
    int k;
    int j;
    long idx;
    double sin;
    double sout;
};

struct IntersectionList {
    struct Intersection *buf;
    int len;
    int cap;
};

struct Grid 
{
    int geom;
    int Nz;
    int Nr;
    int *Np;
    int *idx;
    int Nq;
    double *rjph;
    double *zkph;
    double *piph;
    double *dphi;
    double *prim;
    double t;
    double pmax;
    long Ntot;
};

struct Grid grid_new(char *filename);
void grid_load(struct Grid *g, char filename[]);
void grid_free(struct Grid *g);
void grid_print(struct Grid *g);

void ray_print(struct Ray *r);
void ray_xyz(struct Ray *r, double s, double x[3]);

void grid_pos(struct Grid *g, double xyz[3], double x[3]);
void grid_pos_cyl(struct Grid *g, double xyz[3], double x[3]);

void grid_ray_bound(struct Grid *g, struct Ray *r, double *sin, double *sout);
void grid_ray_bound_cyl(struct Grid *g, struct Ray *r,
                        double *sin, double *sout);

long grid_find_cell_at_pos(struct Grid *g, double xyz[3], int *k, int *j);

struct IntersectionList intersectList_new(int cap);
void intersectList_resize(struct IntersectionList *l, int new_cap);
void intersectList_push(struct IntersectionList *l, struct Intersection *i);
void intersectList_free(struct IntersectionList *l);

struct IntersectionList intersect(struct Grid *g, struct Ray *r);
void refreshIntersections(struct Grid *g, struct Ray *rays,
                            struct IntersectionList *l, int N);
void render(struct Grid *g, struct Ray *rays, struct IntersectionList *ray_hits,
            float *image, int Nx, int Ny,
            void (*transferFunc)(double *, double *, int, float *));

void safe_free(void **p);
int search_sorted(double *arr, int N, double x);
void cross(double a[3], double b[3], double c[3]);
void normalize(double v[3]);
void look_at(double o[3], double target[3], double n[3]);

struct Ray * ray_generate_persp(int Nx, int Ny, double *cam_xyz, double *cam_nz,
                                double *cam_nx, double fov_y_deg);
struct Ray * ray_generate_ortho(int Nx, int Ny, double *cam_xyz, double *cam_nz,
                                double *cam_nx, double dy);

#endif
