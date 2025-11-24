
#include "../disco.h"

static double y_max = 0.0;

void setGeometryParams( struct domain * theDomain ){
    y_max = theDomain->phi_max;
}

double get_dp( double yp , double ym ){
    double dy = yp-ym;
    while( dy < 0.0 ) dy += y_max;
    while( dy > y_max) dy -= y_max;
    return(dy);
}

double get_signed_dp( double yp , double ym ){
    double dy = yp-ym;
    while( dy <-.5*y_max ) dy += y_max;
    while( dy > .5*y_max ) dy -= y_max;
    return(dy);
}

double get_centroid(double xp, double xm, int dim)
{
    return 0.5*(xp+xm);
}

double get_dL( const double * xp , const double * xm , int dim ){
    if( dim==0 )
    {
        double dy = get_dp(xp[1], xm[1]);
        return dy;
    }
    else if( dim==1 )
        return xp[0]-xm[0];
    else 
        return xp[2]-xm[2];
}

double get_dA( const double * xp , const double * xm , int dim ){
    double dx = xp[0]-xm[0];
    double dy = get_dp(xp[1], xm[1]);
    double dz = xp[2]-xm[2];
    if( dim==0 )
        return dx*dz;
    else if( dim==1 )
        return dy*dz;
    else
        return dx*dy;
}

double get_dV( const double * xp , const double * xm ){
    double dx = xp[0]-xm[0];
    double dy = get_dp(xp[1], xm[1]);
    double dz = xp[2]-xm[2];
    return dx*dy*dz;
}

double get_scale_factor( const double * x, int dim)
{
    return 1.0;
}

double get_vol_element(const double *x)
{
    return 1.0;
}

void get_xyz(const double *x, double *xyz)
{
    xyz[0] = x[0];
    xyz[1] = x[1];
    xyz[2] = x[2];
}

void get_rpz(const double *x, double *rpz)
{
    rpz[0] = sqrt(x[0]*x[0]+x[1]*x[1]);
    rpz[1] = atan2(x[1],x[0]);
    rpz[2] = x[2];
}

void get_coords_from_xyz(const double *xyz, double *x)
{
    x[0] = xyz[0];
    x[1] = xyz[1];
    x[2] = xyz[2];
}

void get_coords_from_rpz(const double *rpz, double *x)
{
    x[0] = rpz[0] * cos(rpz[1]);
    x[1] = rpz[0] * sin(rpz[1]);
    x[2] = rpz[2];
}

void get_vec_rpz(const double *x, const double *v, double *vrpz)
{
    double r = sqrt(x[0]*x[0]+x[1]*x[1]);
    double cp = x[0]/r;
    double sp = x[1]/r;
    vrpz[0] =  cp*v[0] + sp*v[1];
    vrpz[1] = -sp*v[0] + cp*v[1];
    vrpz[2] = v[2];
}

void get_vec_from_rpz(const double *x, const double *vrpz, double *v)
{
    double r = sqrt(x[0]*x[0]+x[1]*x[1]);
    double cp = x[0]/r;
    double sp = x[1]/r;
    v[0] = cp*vrpz[0] - sp*vrpz[1];
    v[1] = sp*vrpz[0] + cp*vrpz[1];
    v[2] = vrpz[2];
}

void get_vec_xyz(const double *x, const double *v, double *vxyz)
{
    vxyz[0] = v[0];
    vxyz[1] = v[1];
    vxyz[2] = v[2];
}
void get_vec_from_xyz(const double *x, const double *vxyz, double *v)
{
    v[0] = vxyz[0];
    v[1] = vxyz[1];
    v[2] = vxyz[2];
}

void geom_grad(const double *prim, double *grad, const double *xp, const double *xm, 
                double PLM, int dim, int LR)
{
    printf("Geometric gradient called on non-geometric boundary.\n");
    printf("--Cartesian setup has NO geometric boundaries.\n");
    return;
}

void geom_polar_vec_adjust(const double *xp, const double *xm, double *fac)
{
    fac[0] = 1.0;
    fac[1] = 1.0;
    fac[2] = 1.0;
}

