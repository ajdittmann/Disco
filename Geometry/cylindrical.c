
#include "../paul.h"

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
    if(dim == 1)
        return 2.0*(xp*xp + xp*xm + xm*xm) / (3.0*(xp+xm));
    else
        return 0.5*(xp+xm);
}

double get_dL( const double * xp , const double * xm , int dim ){
    double r = .5*(xp[0]+xm[0]);
    double dphi = get_dp(xp[1], xm[1]);
    if(dim == 0)
        return r*dphi;
    else if(dim == 1)
        return xp[0]-xm[0];
    else
        return xp[2]-xm[2];
}

double get_dA( const double * xp , const double * xm , int dim ){
    double r  = .5*(xp[0]+xm[0]);
    double dr   = xp[0]-xm[0];
    double dphi = get_dp(xp[1], xm[1]);
    double dz   = xp[2]-xm[2];

    if(dim == 0)
        return dr*dz;
    else if(dim == 1)
        return r*dphi*dz;
    else
        return r*dr*dphi;
}

double get_dV( const double * xp , const double * xm ){
    double r  = .5*(xp[0]+xm[0]);
    double dr   = xp[0]-xm[0];
    double dphi = get_dp(xp[1],xm[1]);
    double dz   = xp[2]-xm[2];

    return( r*dr*dphi*dz );
}

double get_scale_factor( const double * x, int dim)
{
    if(dim == 0)
        return x[0];
    return 1.0;
}

double get_vol_element(const double *x)
{
    return x[0];
}

void get_xyz(const double *x, double *xyz)
{
    xyz[0] = x[0] * cos(x[1]);
    xyz[1] = x[0] * sin(x[1]);
    xyz[2] = x[2];
}

void get_rpz(const double *x, double *rpz)
{
    rpz[0] = x[0];
    rpz[1] = x[1];
    rpz[2] = x[2];
}

void get_coords_from_xyz(const double *xyz, double *x)
{
    x[0] = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
    x[1] = atan2(xyz[1],xyz[0]);
    x[2] = xyz[2];
}

void get_coords_from_rpz(const double *rpz, double *x)
{
    x[0] = rpz[0];
    x[1] = rpz[1];
    x[2] = rpz[2];
}

void get_vec_rpz(const double *x, const double *v, double *vrpz)
{
    UNUSED(x);

    vrpz[0] = v[0];
    vrpz[1] = v[1];
    vrpz[2] = v[2];
}

void get_vec_from_rpz(const double *x, const double *vrpz, double *v)
{
    UNUSED(x);

    v[0] = vrpz[0];
    v[1] = vrpz[1];
    v[2] = vrpz[2];
}

void get_vec_xyz(const double *x, const double *v, double *vxyz)
{
    double phi = x[1];
    double cp = cos(phi);
    double sp = sin(phi);

    vxyz[0] = cp*v[0] - sp*v[1];
    vxyz[1] = sp*v[0] + cp*v[1];
    vxyz[2] = v[2];
}

void get_vec_from_xyz(const double *x, const double *vxyz, double *v)
{
    double phi = x[1];
    double cp = cos(phi);
    double sp = sin(phi);

    v[0] =  cp*vxyz[0] + sp*vxyz[1];
    v[1] = -sp*vxyz[0] + cp*vxyz[1];
    v[2] = vxyz[2];
}

void geom_grad(const double *prim, double *grad, const double *xp, const double *xm, 
                double PLM, int dim, int LR)
{
    UNUSED(PLM);

    if(dim !=1 || LR != 0)
    {
        printf("Geometric gradient called on non-geometric boundary\n");
        printf("--Cylindrical setup only has geometric boundary at r=0.\n");
        return;
    }
    if(xp[0] < 0.0 || fabs(xm[0]) > 1.0e-10*xp[0])
    {
        printf("Geometric gradient called on cell with rm = %.le (!= 0)\n",
                xm[0]);
        return;
    }

    int q;
    double r = get_centroid(xp[0], xm[0], 1);
    for(q = 0; q<NUM_Q; q++)
    {
        /*
        if(q == URR || (NUM_C>BZZ && q==BRR))
        {
            double SL = prim[q]/r;
            double S = grad[q];
            if( S*SL < 0.0 )
                grad[q] = 0.0; 
            else if( fabs(PLM*SL) < fabs(S) )
                grad[q] = PLM*SL;
        }
        */
        if(q == UPP)
            grad[q] = -prim[UPP] / r;
        else
            grad[q] = 0.0;
    }
}

void geom_polar_vec_adjust(const double *xp, const double *xm, double *fac)
{
    double dphi = get_dp(xp[1], xm[1]);
    double adjust = sin(0.5*dphi) / (0.5*dphi);

    fac[0] = adjust;
    fac[1] = adjust;
    fac[2] = 1.0;
}

#if ENABLE_CART_INTERP
void geom_interpolate(const double *prim, const double *gradr,
                      const double *gradp, const double *gradz,
                      const double *x, double dr, double dphi, double dz,
                      double * primI, double w)
{

    int q;
    for(q=0; q<NUM_Q; q++)
        primI[q] = prim[q] + dphi * gradp[q] + dr*gradr[q] + dz*gradz[q];

    double r = x[0];
    double sdp = sin(dphi);
    double cdp = cos(dphi);

    // Amend velocity interpolation to include Christoffel terms
    // so interpolation uses covariant derivative instead of coordinate

    // the dr and dphi computed from the cartesian dx and dy
    double dr_cart = dr * cdp - r * (1.0 - cdp);
    double dphi_cart = (r + dr) * sdp / r;

    //Covariant derivatives
    double DurDr = gradr[URR];
    double DurDp = gradp[URR] - r*prim[UPP];
    double DupDr = gradr[UPP] + prim[UPP] / r;
    double DupDp = gradp[UPP] + prim[URR] / r;

    //Extrapolate with covariant derivatives
    double ur_pt = prim[URR] + dr_cart*DurDr + dphi_cart*DurDp + dz*gradz[URR];
    double up_pt = prim[UPP] + dr_cart*DupDr + dphi_cart*DupDp + dz*gradz[UPP];


    primI[URR] = (1.0-w) * primI[URR]
                    + w * ( cdp * ur_pt + r*sdp * up_pt);
    primI[UPP] = (1.0-w) * primI[UPP]
                    + w * (-sdp * ur_pt + r*cdp * up_pt) / (r + dr);
}

void geom_rebase_to_cart(const double *prim, const double *x, double *cartPrim)
{
    int q;
    for(q=0; q<NUM_Q; q++)
        cartPrim[q] = prim[q];

    double r = x[0];
    double cp = cos(x[1]);
    double sp = sin(x[1]);
    cartPrim[URR] = cp * prim[URR] - r * sp * prim[UPP];
    cartPrim[UPP] = sp * prim[URR] + r * cp * prim[UPP];
}

void geom_rebase_from_cart(const double *cartPrim, const double *x,
                           double *prim)
{
    int q;
    for(q=0; q<NUM_Q; q++)
        prim[q] = cartPrim[q];

    double r = x[0];
    double cp = cos(x[1]);
    double sp = sin(x[1]);
    prim[URR] =  cp * cartPrim[URR] + sp * cartPrim[UPP];
    prim[UPP] = (-sp * cartPrim[URR] + cp * cartPrim[UPP]) / r;
}

void geom_gradCart_to_grad(const double *cartGrad, const double *prim,
                           const double *x, double *grad, int dim)
{
    int q;
    for(q=0; q<NUM_Q; q++)
        grad[q] = cartGrad[q];

    double r = x[0];
    double cp = cos(x[1]);
    double sp = sin(x[1]);
    // grad[Uxx] is now the covariant derivative of Uxx in polar coords
    grad[URR] =  cp * cartGrad[URR] + sp * cartGrad[UPP];
    grad[UPP] = (-sp * cartGrad[URR] + cp * cartGrad[UPP]) / r;

    if(dim == 0) // phi
    {
        grad[URR] -= -r * prim[UPP];  // G^r_pp * v^p
        grad[UPP] -= prim[URR] / r;   // G^p_pr * v^r
    }
    else if (dim == 1)
    {
        grad[UPP] -= prim[UPP] / r;   // G^p_rp * v^p
    }
}

void geom_cart_interp_grad_trans(const double *primL, const double *primR,
                                 const double *gradpL, const double *gradpR,
                                 double dpL, double dpR, const double *x,
                                 double dxL, double dxR, 
                                 double *gradCIL, double *gradCIR,
                                 int dim)
{
    /*
     * In this function, dxL & dxR are SIGNED.  That is, dxL < 0.0
     * and dxR > 0.0
     */

    int q;

    double idx = 1.0 / (dxR - dxL);
    for(q=0; q<NUM_Q; q++)
    {
        gradCIL[q] = idx * (primR[q] + dpR * gradpR[q]
                            - (primL[q] + dpL * gradpL[q]));
        gradCIR[q] = gradCIL[q];
    }

    double xL[3] = {x[0], x[1]-dpL, x[2]};
    double xR[3] = {x[0], x[1]-dpR, x[2]};
    if(dim == 1)
    {
        xL[0] += dxL;
        xR[0] += dxR;
    }
    else if(dim == 2)
    {
        xL[2] += dxL;
        xR[2] += dxR;
    }

    // Defining a local cartesian frame centered on the face
    // In this frame the face is on the x-axis

    double rL = xL[0];
    double rR = xR[0];
    double sL = sin(-dpL);  //sin(xL[1]);
    double cL = cos(-dpL);  //cos(xL[1]);
    double sR = sin(-dpR);  //sin(xR[1]);
    double cR = cos(-dpR);  //cos(xR[1]);

    double uxL = cL * primL[URR] - rL * sL * primL[UPP];
    double uyL = sL * primL[URR] + rL * cL * primL[UPP];
    double uxR = cR * primR[URR] - rR * sR * primR[UPP];
    double uyR = sR * primR[URR] + rR * cR * primR[UPP];

    double DurDpL = gradpL[URR] - rL * primL[UPP];
    double DupDpL = gradpL[UPP] + primL[URR] / rL;
    double DurDpR = gradpR[URR] - rR * primR[UPP];
    double DupDpR = gradpR[UPP] + primR[URR] / rR;

    double duxdpL = cL * DurDpL - rL * sL * DupDpL;
    double duydpL = sL * DurDpL + rL * cL * DupDpL;
    double duxdpR = cR * DurDpR - rR * sR * DupDpR;
    double duydpR = sR * DurDpR + rR * cR * DupDpR;

    double dpR_cart = sin(dpR);
    double dpL_cart = sin(dpL);
    double drR_cart = rR * (cos(dpR) - 1);
    double drL_cart = rL * (cos(dpL) - 1);

    dpR_cart += drR_cart * sR / (rR * cR);
    dpL_cart += drL_cart * sL / (rL * cL);

    double gradUX = idx * (uxR + dpR_cart * duxdpR
                            - (uxL + dpL_cart * duxdpL));
    double gradUY = idx * (uyR + dpR_cart * duydpR
                            - (uyL + dpL_cart * duydpL));

    gradUX /= (1 - idx*(drR_cart/cR-drL_cart/cL));
    gradUY /= (1 - idx*(drR_cart/cR-drL_cart/cL));

    double gradUXL = (gradUX + duxdpL * sL / rL) / cL;
    double gradUYL = (gradUY + duydpL * sL / rL) / cL;
    double gradUXR = (gradUX + duxdpR * sR / rR) / cR;
    double gradUYR = (gradUY + duydpR * sR / rR) / cR;

    gradCIL[URR] =  cL * gradUXL + sL * gradUYL;
    gradCIL[UPP] = (-sL * gradUXL + cL * gradUYL) / rL;
    gradCIR[URR] =  cR * gradUXR + sR * gradUYR;
    gradCIR[UPP] = (-sR * gradUXR + cR * gradUYR) / rR;

    //
    // gradCI are now the *covariant* derivatives at L and R.
    // These are actually what we want, as the gradients are initialized with
    // the 0th order values, which are just the negative of the
    // christoffel-symbols part.  When gradCI are added, then, 
    // we'll be left with the desired coordinate derivatives.
    //

    /*
    if(dim == 1)
    {
        gradCIL[UPP] -= primL[UPP] / rL;
        gradCIR[UPP] -= primR[UPP] / rR;
    }
    */
}
#endif
