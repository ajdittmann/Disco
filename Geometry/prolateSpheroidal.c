
#include "../disco.h"

/*
 * Prolate Spheroidal Coordinates mu, phi, nu
 *
 * a = 1
 *
 * ds^2 = a^2[(sinh^2(mu) + sin^2(nu))(dmu^2 + dnu^2)
 *          + sinh^2(mu) sin^2(nu) * dphi^2
 *
 * mu << a ==>  mu ~ r/a, nu ~ z/a
 * mu >> a ==>  mu ~ R/a, nu ~ theta
 */

static double phi_max = 0.0;

void setGeometryParams( struct domain * theDomain ){
   phi_max = theDomain->phi_max;
}

double get_dp( double phip , double phim ){
   double dp = phip-phim;
   while( dp < 0.0 ) dp += phi_max;
   while( dp > phi_max) dp -= phi_max;
   return(dp);
}

double get_signed_dp( double phip , double phim ){
   double dp = phip-phim;
   while( dp <-.5*phi_max ) dp += phi_max;
   while( dp > .5*phi_max ) dp -= phi_max;
   return(dp);
}

double get_centroid(double xp, double xm, int dim)
{
    return 0.5*(xp+xm);
}

double get_dL( const double * xp , const double * xm , int dim ){
    double mu = 0.5*(xp[0]+xm[0]);
    double nu = 0.5*(xp[2]+xm[2]);
    double s_mu = sinh(mu);
    double s_nu = sin(nu);
    if(dim == 0)
    {
        double dphi = get_dp(xp[1], xm[1]);
        return s_mu*s_nu*dphi;
    }
    else if(dim == 1)
        return sqrt(s_mu*s_mu + s_nu*s_nu) * (xp[0] - xm[0]);
    else
        return sqrt(s_mu*s_mu + s_nu*s_nu) * (xp[2] - xm[2]);
}

double get_dA( const double * xp , const double * xm , int dim ){
    
    double mu = 0.5*(xp[0]+xm[0]);
    double nu = 0.5*(xp[2]+xm[2]);
    double s_mu = sinh(mu);
    double s_nu = sin(nu);

    if(dim == 0)
    {
        double dmu = xp[0]-xm[0];
        double dnu = xp[2]-xm[2];
        return (s_mu*s_mu + s_nu*s_nu)*dmu*dnu;
    }
    else if(dim == 1)
    {
        double dnu = xp[2]-xm[2];
        double dphi = get_dp(xp[1], xm[1]);
        return s_mu*s_nu*sqrt(s_mu*s_mu + s_nu*s_nu)*dnu*dphi;
    }
    else
    {
        double dmu = xp[0]-xm[0];
        double dphi = get_dp(xp[1], xm[1]);
        return s_mu*s_nu*sqrt(s_mu*s_mu + s_nu*s_nu)*dmu*dphi;
    }
}

double get_dV( const double * xp , const double * xm ){
    double mu = 0.5*(xp[0]+xm[0]);
    double nu = 0.5*(xp[2]+xm[2]);
    double s_mu = sinh(mu);
    double s_nu = sin(nu);
    double dmu = xp[0]-xm[0];
    double dnu = xp[2]-xm[2];
    double dphi = get_dp(xp[1], xm[1]);

    return (s_mu*s_mu + s_nu*s_nu) * s_mu*s_nu*dmu*dnu*dphi;
}

double get_scale_factor( const double * x, int dim)
{
    double s_mu = sinh(x[0]);
    double s_nu = sin(x[2]);
    if(dim == 0)
        return s_mu*s_nu;
    else
        return sqrt(s_mu*s_mu + s_nu*s_nu);
}

double get_vol_element(const double *x)
{
    double s_mu = sinh(x[0]);
    double s_nu = sin(x[2]);
    return (s_mu*s_mu + s_nu*s_nu) * s_mu*s_nu;
}

void get_xyz(const double *x, double *xyz)
{
    double s_mu = sinh(x[0]);
    double s_nu = sin(x[2]);
    xyz[0] = s_mu * s_nu * cos(x[1]);
    xyz[1] = s_mu * s_nu * sin(x[1]);
    xyz[2] = cosh(x[0]) * cos(x[2]);
}

void get_rpz(const double *x, double *rpz)
{
    rpz[0] = sinh(x[0]) * sin(x[2]);
    rpz[1] = x[1];
    rpz[2] = cosh(x[0]) * cos(x[2]);
}

void get_coords_from_xyz(const double *xyz, double *x)
{
    double r = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
    double z = xyz[2];
    double Rp = sqrt(r*r + (z-1)*(z-1));
    double Rm = sqrt(r*r + (z+1)*(z+1));

    x[0] = acosh(0.5*(Rm + Rp));
    x[1] = atan2(xyz[1],xyz[0]);
    x[2] = acos(0.5*(Rm - Rp));
}

void get_coords_from_rpz(const double *rpz, double *x)
{
    double r = rpz[0];
    double z = rpz[2];
    double Rp = sqrt(r*r + (z-1)*(z-1));
    double Rm = sqrt(r*r + (z+1)*(z+1));

    x[0] = acosh(0.5*(Rm + Rp));
    x[1] = rpz[1];
    x[2] = acos(0.5*(Rm - Rp));
}

void get_vec_rpz(const double *x, const double *v, double *vrpz)
{
    double s_mu = sinh(x[0]);
    double c_mu = cosh(x[0]);
    double s_nu = sin(x[2]);
    double c_nu = cos(x[2]);
    double h = sqrt(s_mu*s_mu + s_nu*s_nu);

    vrpz[0] =  (c_mu*s_nu*v[0] + s_mu*c_nu*v[2]) / h;
    vrpz[1] = v[1];
    vrpz[2] =  (s_mu*c_nu*v[0] - c_mu*s_nu*v[2]) / h;
}

void get_vec_from_rpz(const double *x, const double *vrpz, double *v)
{
    double s_mu = sinh(x[0]);
    double c_mu = cosh(x[0]);
    double s_nu = sin(x[2]);
    double c_nu = cos(x[2]);
    double h = sqrt(s_mu*s_mu + s_nu*s_nu);

    v[0] = (c_mu*s_nu * vrpz[0] + s_mu*c_nu * vrpz[2]) / h;
    v[1] = vrpz[1];
    v[2] = (s_mu*c_nu * vrpz[0] - c_mu*s_nu * vrpz[2]) / h;
}

void get_vec_xyz(const double *x, const double *v, double *vxyz)
{
    double s_mu = sinh(x[0]);
    double c_mu = cosh(x[0]);
    double s_nu = sin(x[2]);
    double c_nu = cos(x[2]);
    double s_p = sin(x[1]);
    double c_p = cos(x[1]);
    double h = sqrt(s_mu*s_mu + s_nu*s_nu);

    double vr = (c_mu*s_nu*v[0] + s_mu*c_nu*v[2]) / h;
    double vp = v[1];
    vxyz[0] = c_p * vr - s_p * vp;
    vxyz[1] = s_p * vr + c_p * vp;
    vxyz[2] =  (s_mu*c_nu*v[0] - c_mu*s_nu*v[2]) / h;
}

void get_vec_from_xyz(const double *x, const double *vxyz, double *v)
{
    double s_mu = sinh(x[0]);
    double c_mu = cosh(x[0]);
    double s_nu = sin(x[2]);
    double c_nu = cos(x[2]);
    double s_p = sin(x[1]);
    double c_p = cos(x[1]);
    double h = sqrt(s_mu*s_mu + s_nu*s_nu);

    double vr =  c_p * vxyz[0] + s_p * vxyz[1];
    double vp = -s_p * vxyz[0] + c_p * vxyz[1];

    v[0] = (c_mu*s_nu * vr + s_mu*c_nu * vxyz[2]) / h;
    v[1] = vp;
    v[2] = (s_mu*c_nu * vr - c_mu*s_nu * vxyz[2]) / h;
}

void geom_grad(const double *prim, double *grad, const double *xp, const double *xm, 
                double PLM, int dim, int LR)
{
    if(!(dim==2 || (dim==1 && LR ==0)))
    {
        printf("Geometric gradient called on non-geometric boundary\n");
        printf("--Prolate Spheroidal setup only has geometric boundary at mu=0"
               " and nu=0,pi.\n");
        return;
    }
    if(dim == 1)
    {
        if(xp[0] < 0.0 || fabs(xm[0]) > 1.0e-10*xp[0])
        {
            printf("Geometric gradient called on cell with rm = %le (!= 0)\n",
                    xm[0]);
            return;
        }

        int q;
        double r = get_centroid(xp[0], xm[0], 1);
        for(q = 0; q<NUM_Q; q++)
        {
            if(q == URR || (NUM_C>BZZ && q==BRR))
            {
                double SL = prim[q]/r;
                double S = grad[q];
                if( S*SL < 0.0 )
                    grad[q] = 0.0; 
                else if( fabs(PLM*SL) < fabs(S) )
                    grad[q] = PLM*SL;
            }
            else
                grad[q] = 0.0;
        }
    }
    if(dim == 2)
    {
        if(LR == 0 && (xm[2] < 0.0 || fabs(xm[2]) > 1.0e-10*xp[2]))
        {
            printf("Geometric gradient called on cell with"
                    " theta_m = %le (!= 0)\n", xm[2]);
            return;
        }
        else if(LR == 1 && (fabs(xp[2]-M_PI) > 1.0e-10*(xp[2]-xm[2])))
        {
            printf("Geometric gradient called on cell with"
                    " theta_p = %le (!= pi)\n", xp[2]);
            return;
        }

        int q;
        double th = get_centroid(xp[2], xm[2], 2);
        double dth = LR==0 ? -th : th-M_PI; //dth is negative on th=pi to get
                                            // the correct sign on SL.
        for(q = 0; q<NUM_Q; q++)
        {
            if(q == UZZ || (NUM_C>BZZ && q==BZZ))
            {
                double SL = prim[q]/dth;
                double S = grad[q];
                if( S*SL < 0.0 )
                    grad[q] = 0.0; 
                else if( fabs(PLM*SL) < fabs(S) )
                    grad[q] = PLM*SL;
            }
            else
                grad[q] = 0.0;
        }
    }
}

void geom_polar_vec_adjust(const double *xp, const double *xm, double *fac)
{
    //TODO: IMPLEMENT THIS!!!
    
    fac[0] = 1.0;
    fac[1] = 1.0;
    fac[2] = 1.0;
}
