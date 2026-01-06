
#include "../disco.h"

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
    if(xp == xm)
        return xp;

    if(dim == 1)
    {
        double x2_3 = xp*xp + xp*xm + xm*xm;
        return 3.0*(xp*x2_3 + xm*xm*xm) / (4.0*x2_3);
    }
    else if(dim == 2)
    {
        double cp = cos(xp);
        double sp = sin(xp);
        double cm = cos(xm);
        double sm = sin(xm);
        return (xm*cm-xp*cp+sp-sm) / (cm-cp);
    }
    else
        return 0.5*(xp+xm);
}

double get_dL( const double * xp , const double * xm , int dim ){
    double r = .5*(xp[0]+xm[0]);
    if(dim == 0)
    {
        double th = .5*(xp[2]+xm[2]);
        double sth = sin(th);
        double dphi = get_dp(xp[1], xm[1]);
        return r*sth*dphi;
    }
    else if(dim == 1)
        return xp[0]-xm[0];
    else
        return r*(xp[2]-xm[2]);
}

double get_dA( const double * xp , const double * xm , int dim ){
    
    double r  = .5*(xp[0]+xm[0]);

    if(dim == 0)
    {
        double dr = xp[0]-xm[0];
        double dth = xp[2]-xm[2];
        return r*dr*dth;
    }
    else if(dim == 1)
    {
        double sinth = sin(0.5*(xp[2]+xm[2]));
        double sindth = 2*sin(0.5*(xp[2]-xm[2]));
        double dphi = get_dp(xp[1], xm[1]);
        return r*r*sinth*sindth*dphi;
    }
    else
    {
        double dr = xp[0]-xm[0];
        double sinth = sin(0.5*(xp[2]+xm[2]));
        double dphi = get_dp(xp[1], xm[1]);
        return r*sinth*dr*dphi;
    }
}

double get_dV( const double * xp , const double * xm ){
    double r2 = (xp[0]*xp[0]+xp[0]*xm[0]+xm[0]*xm[0])/3.0;
    double dr   = xp[0]-xm[0];
    double dphi = get_dp(xp[1],xm[1]);
    double sinth = sin(0.5*(xp[2]+xm[2]));
    double sindth = 2*sin(0.5*(xp[2]-xm[2]));

    return r2*sinth*dr*sindth*dphi;
}

double get_scale_factor( const double * x, int dim)
{
    if(dim == 0)
        return x[0]*sin(x[2]);
    else if(dim == 2)
        return x[0];
    return 1.0;
}

double get_vol_element(const double *x)
{
    return x[0]*x[0]*sin(x[2]);
}

void get_xyz(const double *x, double *xyz)
{
    double sinth = sin(x[2]);
    xyz[0] = x[0] * cos(x[1]) * sinth ;
    xyz[1] = x[0] * sin(x[1]) * sinth;
    xyz[2] = x[0] * cos(x[2]);
}

void get_rpz(const double *x, double *rpz)
{
    rpz[0] = x[0]*sin(x[2]);
    rpz[1] = x[1];
    rpz[2] = x[0]*cos(x[2]);
}

void get_coords_from_xyz(const double *xyz, double *x)
{
    double r = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
    x[0] = r;
    x[1] = atan2(xyz[1],xyz[0]);
    x[2] = acos(xyz[2] / r);
}

void get_coords_from_rpz(const double *rpz, double *x)
{
    double r = sqrt(rpz[0]*rpz[0]+rpz[2]*rpz[2]);
    x[0] = r;
    x[1] = rpz[1];
    x[2] = acos(rpz[2] / r);
}

void get_vec_rpz(const double *x, const double *v, double *vrpz)
{
    double st = sin(x[2]);
    double ct = cos(x[2]);
    vrpz[0] =  ct*v[2] + st*v[0];
    vrpz[1] = v[1];
    vrpz[2] = -st*v[2] + ct*v[0];
}

void get_vec_from_rpz(const double *x, const double *vrpz, double *v)
{
    double st = sin(x[2]);
    double ct = cos(x[2]);
    v[0] = ct*vrpz[2] + st*vrpz[0];
    v[1] = vrpz[1];
    v[2] = -st*vrpz[2] + ct*vrpz[0];
}

void get_vec_xyz(const double *x, const double *v, double *vxyz)
{
    double cp = cos(x[1]);
    double sp = sin(x[1]);
    double ct = cos(x[2]);
    double st = sin(x[2]);
    
    //Taken straight from the back cover of Griffith's E&M
    vxyz[0] = st*cp*v[0] + ct*cp*v[2] - sp*v[1];
    vxyz[1] = st*sp*v[0] + ct*sp*v[2] + cp*v[1];
    vxyz[2] = ct*v[0] - st*v[2];
}

void get_vec_from_xyz(const double *x, const double *vxyz, double *v)
{
    double cp = cos(x[1]);
    double sp = sin(x[1]);
    double ct = cos(x[2]);
    double st = sin(x[2]);
    
    //Taken straight from the back cover of Griffith's E&M

    v[0] =  st*cp*vxyz[0] + st*sp*vxyz[1] + ct*vxyz[2];
    v[1] = -sp*vxyz[0] + cp*vxyz[1];
    v[2] = ct*cp*vxyz[0] + ct*sp*vxyz[1] - st*vxyz[2];
}

void geom_grad(const double *prim, double *grad, const double *xp, const double *xm, 
                double PLM, int dim, int LR)
{
    if(!(dim==2 || (dim==1 && LR ==0)))
    {
        printf("Geometric gradient called on non-geometric boundary\n");
        printf("--Spherical setup only has geometric boundary at r=0 and th=0,pi.\n");
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
            printf("Geometric gradient called on cell with theta_m = %le (!= 0)\n",
                    xm[2]);
            return;
        }
        else if(LR == 1 && (fabs(xp[2]-M_PI) > 1.0e-10*(xp[2]-xm[2])))
        {
            printf("Geometric gradient called on cell with theta_p = %le (!= pi)\n",
                    xp[2]);
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
