#include "paul.h"
#include "geometry.h"
#include "planet.h"

static int grav2D = 0;
static int polar_sources_r = 0;
static int polar_sources_p = 0;
static int polar_sources_z = 0;

void setGravParams( struct domain * theDomain ){

   grav2D = theDomain->theParList.grav2D;
   if(strcmp(GEOMETRY, "cylindrical") == 0
             && theDomain->theParList.NoBC_Rmin == 1)
   {
       polar_sources_r = 1;
       polar_sources_p = 1;
   }
   if(strcmp(GEOMETRY, "spherical") == 0)
   {
       if(theDomain->theParList.NoBC_Rmin == 1)
       {
           polar_sources_r = 1;
           polar_sources_p = 1;
       }
       if(theDomain->theParList.NoBC_Zmin == 1
          || theDomain->theParList.NoBC_Zmax == 1)
       {
           polar_sources_p = 1;
           polar_sources_z = 1;
       }
   }

}

double phigrav( double M , double r , double eps , int type)
{
    if(type == PLPOINTMASS)
    {
        return( M / sqrt( r*r + eps*eps) ) ;
    }
    else if(type == PLPW)
    {
        return M / (r - 2*M);
    }
    else if(type == PLWEGGC)
    {
        //potential "C" from Wegg 2012, ApJ 749
        double sq6 = sqrt(6);
        double Alpha = -4*(2.0+sq6)/3;
        double Rx = M*(4.0*sq6 - 9);
        double Ry = -4*M*(2*sq6 - 3.0)/3.0;
        return -Alpha*M/r - (1-Alpha)*M/(r-Rx) - M*Ry/(r*r);
    }
    else if(type == PLSURFACEGRAV)
    {
        return M*r; // M is gravitational acceleration
                    // only makes sense if grav2D is on
    }
    else if(type == PLQUAD)
    {
        return 0.5*M*r*r;
    }
    else if(type == PLSPLINE)
    {
        eps = eps*2.8;
        double u = r/eps;
        double u2 = u*u;
        double u3 = u2*u;
        double u4 = u2*u2;
        double u5 = u2*u3;
        double val;
        if (u<0.5)
            val = 16.*u2/3. - 48.*u4/5. + 32.*u5/5. - 14./5.;
        else if (u < 1.0)
            val = 1./(15.*u) + 32.*u2/3. - 16.*u3 + 48.*u4/5. - 32.*u5/15. - 3.2;
        else
            val = -1./u ;

        return -1*M*val/eps;
    }
    else if(type == PLUNIFORM)
    {
        return 0.0;
    }
    else if(type == PLLINEARX)
    {
        return 0.0;
    }
    return 0.0;
}

double fgrav( double M , double r , double eps , int type)
{
    if(type == PLPOINTMASS)
    {
        double R = sqrt(r*r + eps*eps);
        return( M*r / (R*R*R) );
    }
    else if(type == PLPW)
    {
        return M / ((r-2*M)*(r-2*M));
    }
    else if(type == PLSURFACEGRAV)
    {
        return M; // M is gravitational acceleration
                  // only makes sense if grav2D is on
    }
    else if(type == PLQUAD)
    {
        return M * r;
    }
    else if(type == PLWEGGC)
    {
        //potential "C" from Wegg 2012, ApJ 749
        double sq6 = sqrt(6);
        double Alpha = -4*(2.0+sq6)/3;
        double Rx = M*(4.0*sq6 - 9);
        double Ry = -4*M*(2*sq6 - 3.0)/3.0;
        return -Alpha*M/(r*r) - (1-Alpha)*M/((r-Rx)*(r-Rx)) - 2*M*Ry/(r*r*r);
    }
    else if(type == PLSPLINE)
    {
        eps = eps*2.8;
        double u = r/eps;
        double u2 = u*u;
        double u3 = u2*u;
        double u4 = u2*u2;
        double val;
        if (u<0.5)
            val = 32.*u/3. - 192.*u3/5. + 32.*u4;
        else if (u < 1.0)
            val = -1./(15.*u2) + 64.*u/3. - 48.*u2 + 192.*u3/5. - 32*u4/3.;
        else
            val = 1./u2;
        return M*val/(eps*eps);
    }
    else if(type == PLUNIFORM)
    {
        return M;
    }
    else if(type == PLLINEARX)
    {
        return M*r;
    }
    return 0.0;
    
}

void adjust_gas( struct planet * pl , double * x , double * prim , double gam ){

   double r   = x[0];
   double phi = x[1];

   double rp = pl->r;
   double pp = pl->phi;
   double cosp = cos(phi);
   double sinp = sin(phi);
   double dx = r*cosp-rp*cos(pp);
   double dy = r*sinp-rp*sin(pp);
   double script_r = sqrt(dx*dx+dy*dy);

   //double z = M_PI*0.5; //x[2];
   double pot = phigrav( pl->M , script_r , pl->eps , pl->type);

   double c2 = gam*prim[PPP]/prim[RHO];
   double factor = 1. + (gam-1.)*pot/c2;

   prim[RHO] *= factor;
   prim[PPP] *= pow( factor , gam );

}

double planetaryPotential(struct planet *pl, const double *xyz)
{
    /*
     * Returns the gravitational potential (with correct sign) at position
     * (r, phi, z) from planet pl.
     */
    
    double xp = pl->xyz[0];
    double yp = pl->xyz[1];
    double zp = pl->xyz[2];

    double z = xyz[2];

    if(grav2D == 1)
        z = zp;
    else if(grav2D == 2)
    {
        xp = xyz[0];
        yp = xyz[1];
        z = 1.0;
    }

    double dx = xyz[0] - xp;
    double dy = xyz[1] - yp;
    double dz = z - zp;
    double script_r = sqrt(dx*dx + dy*dy + dz*dz);

    return -phigrav(pl->M, script_r, pl->eps, pl->type);
}

void planetaryForce( struct planet * pl , const double *xyz, double *Fxyz)
{
    /*
     * Calculates the specific gravitational force (ie. acceleration) from a
     * planet at location xyz. The force on the planet is the -'ve of this
     * value.
     */
    
    double xp = pl->xyz[0];
    double yp = pl->xyz[1];
    double zp = pl->xyz[2];

    double z = xyz[2];

    if(grav2D == 1)
        z = zp;
    else if(grav2D == 2)
    {
        xp = xyz[0];
        yp = xyz[1];
        z = 1.0;
    }

    if(pl->type == PLUNIFORM)
    {
        xp = xyz[0] + 1.0*cos(pl->phi);
        yp = xyz[1] + 1.0*sin(pl->phi);
        zp = xyz[2];
    }
    else if(pl->type == PLLINEARX)
    {
        xp = xyz[0] + xyz[0]*cos(pl->phi);
        yp = xyz[1] + xyz[0]*sin(pl->phi);
        zp = xyz[2];
    }

    double dx = xyz[0] - xp;
    double dy = xyz[1] - yp;
    double dz = z - zp;
    double script_r = sqrt(dx*dx + dy*dy + dz*dz);

    double f1 = -fgrav( pl->M , script_r , pl->eps , pl->type);

    double ir = 1.0 / script_r;

    Fxyz[0] = f1 * dx * ir;
    Fxyz[1] = f1 * dy * ir;
    Fxyz[2] = f1 * dz * ir;
}

void planet_src( struct planet * pl, const double * prim, double * cons,
                const double * xp, const double * xm, const double *xyz,
                double dV, double dt, double *pl_gas_track )
{
   double rho = prim[RHO];
   double vr  = prim[URR];
   double vz  = prim[UZZ];
   double omega = prim[UPP];
   
   //double rp = xp[0];
   //double rm = xm[0];
   //double r = 0.5*(rp+rm);
  
   double x[3];
   get_centroid_arr(xp, xm, x);
   
   double Fxyz[3], F[3];
   planetaryForce( pl, xyz, Fxyz);
   get_vec_from_xyz(x, Fxyz, F);


   // if we're doing polar stuff, F[] gets modified
   if(polar_sources_r || polar_sources_p || polar_sources_z)
   {
       double adjust[3];
       geom_polar_vec_adjust(xp, xm, adjust);
       if(polar_sources_r)
           F[0] *= adjust[0];
       if(polar_sources_p)
           F[1] *= adjust[1];
       if(polar_sources_z)
           F[2] *= adjust[2];

       // Modified F[], so re-calculate Fxyz
       get_vec_xyz(x, F, Fxyz);
   }

   //Put F in the covariant basis. In cylindrical coords, this multiplies
   //F_phi by 'r'. This is the form needed for the source terms and for
   //dotting against the velocities.
   get_vec_covariant(x, F, F);

   double dVdt = dV*dt;

   cons[SRR] += rho*F[0]*dVdt;
   cons[LLL] += rho*F[1]*dVdt;
   cons[SZZ] += rho*F[2]*dVdt;
   cons[TAU] += rho*( F[0]*vr + F[1]*omega + F[2]*vz )*dVdt;

   double rFp = pl->xyz[0] * Fxyz[1] - pl->xyz[1] * Fxyz[0];

   double Phi = planetaryPotential(pl, xyz);

   pl_gas_track[PL_GRV_PX] -= rho*Fxyz[0]*dVdt;
   pl_gas_track[PL_GRV_PY] -= rho*Fxyz[1]*dVdt;
   pl_gas_track[PL_GRV_PZ] -= rho*Fxyz[2]*dVdt;
   pl_gas_track[PL_GRV_JZ] -= rho*rFp*dVdt;
   pl_gas_track[PL_GRV_EGAS] -= rho*(F[0]*vr + F[1]*omega + F[2]*vz)*dVdt;
   pl->Uf += rho*Phi*dV;
}

void copyPlanetsRK( struct domain * theDomain ){
    int pq;

    double *pl_kin = theDomain->pl_kin;
    double *pl_RK_kin = theDomain->pl_RK_kin;
    double *pl_aux = theDomain->pl_aux;
    double *pl_RK_aux = theDomain->pl_RK_aux;

    for(pq=0; pq < theDomain->Npl * NUM_PL_KIN; pq++)
        pl_RK_kin[pq] = pl_kin[pq];

    for(pq=0; pq < theDomain->Npl * NUM_PL_AUX; pq++)
        pl_RK_aux[pq] = pl_aux[pq];
}

void adjustPlanetsRKkin( struct domain * theDomain , double RK ){
    int pq;

    double *pl_kin = theDomain->pl_kin;
    double *pl_RK_kin = theDomain->pl_RK_kin;

    for(pq=0; pq < theDomain->Npl * NUM_PL_KIN; pq++)
        pl_kin[pq] = (1-RK) * pl_kin[pq] + RK * pl_RK_kin[pq];
}

void adjustPlanetsRKaux( struct domain * theDomain , double RK ){
    int pq;

    double *pl_aux = theDomain->pl_aux;
    double *pl_RK_aux = theDomain->pl_RK_aux;
    
    for(pq=0; pq < theDomain->Npl * NUM_PL_AUX; pq++)
        pl_aux[pq] = (1-RK) * pl_aux[pq] + RK * pl_RK_aux[pq];
}

void planet_init_kin(struct planet *pl, double *pl_kin)
{
    /* Set the planet kin[] parameters consistent with the
     * planet's M, r, phi, etc.
     * Analagous to prim2cons() in Hydro.
     */

    pl_kin[PL_M] = pl->M;
    pl_kin[PL_R] = pl->r;
    pl_kin[PL_PHI] = pl->phi;
    pl_kin[PL_Z] = 0.0;
    pl_kin[PL_PR] = pl->M * pl->vr;
    pl_kin[PL_LL] = pl->M * pl->r * pl->r * pl->omega;
    pl_kin[PL_PZ] = 0.0;
    pl_kin[PL_SZ] = 0.0;
    pl_kin[PL_EINT] = 0.0;
}

void zeroAuxPlanets( struct domain *theDomain)
{
    int pq;
    for(pq=0; pq < NUM_PL_AUX * theDomain->Npl; pq++)
    {
        theDomain->pl_aux[pq] = 0.0;
        theDomain->pl_RK_aux[pq] = 0.0;
    }
}

void setupPlanets(struct domain *theDomain)
{
   initializePlanets( theDomain->thePlanets );
   setPlanetsXYZ(theDomain);
  
   int Npl = theDomain->Npl;
   int p;
   for(p=0; p<Npl; p++)
   {
       planet_init_kin(theDomain->thePlanets + p,
                       theDomain->pl_kin + p*NUM_PL_KIN);
   }

   zeroAuxPlanets(theDomain);
   initializePlanetTracking(theDomain);
}

void setPlanetsXYZ(struct domain *theDomain)
{
    int p;
    for(p=0; p<theDomain->Npl; p++)
    {
        struct planet *pl = theDomain->thePlanets + p;
        pl->xyz[0] = pl->r * cos(pl->phi);
        pl->xyz[1] = pl->r * sin(pl->phi);
        pl->xyz[2] = 0.0;
    }
}

void movePlanetsLive(struct domain *theDomain)
{
    /* 
     * Adjusts planet mass, position, and velocity from kin[]
     * values.  Like cons2prim.
     */

    int p;
    int Npl = theDomain->Npl;
    struct planet *thePlanets = theDomain->thePlanets;

    for(p=0; p<Npl; p++)
    {
        struct planet *pl = thePlanets + p;
        double *pl_kin = theDomain->pl_kin + p*NUM_PL_KIN;

        double r = pl_kin[PL_R];
        double M = pl_kin[PL_M];

        pl->M = M;
        pl->r = r;
        pl->phi = pl_kin[PL_PHI];
        pl->vr = pl_kin[PL_PR] / M;
        pl->omega = pl_kin[PL_LL] / (M * r * r);
    }
}

void initializePlanetTracking(struct domain *theDomain)
{
    int Npl = theDomain->Npl;
    int p, pq;
    
    for(p=0; p<Npl; p++)
        theDomain->thePlanets[p].Uf = 0.0;
    
    for(pq=0; pq<Npl * NUM_PL_INTEGRALS; pq++)
        theDomain->pl_gas_track[pq] = 0.0;

    theDomain->planet_gas_track_synced = 0;
}

void updatePlanetsKinAux(struct domain *theDomain, double dt)
{
    /*
     * If the planets are LIVE, then the gas_track integrals *must*
     * have been reduced (summed) over MPI ranks already for this to work.
     *
     * This is like updating cons from the fluxes.
     */
    int p;
    int Npl = theDomain->Npl;
    struct planet *thePlanets = theDomain->thePlanets;
   
    int live_planet = !planet_motion_analytic();

    for(p=0; p<Npl; p++)
    {
        struct planet *pl = thePlanets + p;

        double cosp = cos(pl->phi);
        double sinp = sin(pl->phi);

        double M = pl->M;
        double r = pl->r;
        double omega = pl->omega;
        
        double x = r * cosp;
        double y = r * sinp;
        double vx = pl->vr * cosp - r * pl->omega * sinp;
        double vy = pl->vr * sinp + r * pl->omega * cosp;
        double vz = 0.0;

        double v2 = vx*vx + vy*vy + vz*vz;


        // Compute the gravitational potential and acceleration due to other
        // planets.

        double Phi_ext = 0.0;
        double gx_ext = 0;
        double gy_ext = 0;
        double gz_ext = 0;

        int p2;

        for(p2=0; p2<Npl; p2++)
        {
            if(p == p2)
                continue;
        
            struct planet *pl2 = theDomain->thePlanets + p2;
            double r2 = pl2->r;
            double cosp2 = cos(pl2->phi);
            double sinp2 = sin(pl2->phi);
            double x2 = r2*cosp2;
            double y2 = r2*sinp2;

            double rsep = sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2));

            Phi_ext += -phigrav(pl2->M, rsep, pl2->eps, pl2->type);
            double g = fgrav(pl2->M, rsep, pl2->eps, pl2->type);

            gx_ext += g * (x2-x) / rsep;
            gy_ext += g * (y2-y) / rsep;
            gz_ext += 0;
        }

        // These "*_pl" values are only set if the planet is live.
        // In this case, the "kin" values track the full evolution of the
        // planet.  Otherwise, "kin" only stores the *additional* accumulated
        // offset from accretion

        double vr_pl = 0;
        double om_pl = 0;
        double vz_pl = 0;
        double Fr_pl = 0;
        double rFp_pl = 0;
        double Fz_pl = 0;
        double Fx_pl = 0;
        double Fy_pl = 0;

        if(live_planet)
        {
            vr_pl = pl->vr;
            om_pl = omega;

            double gr = cosp * gx_ext + sinp * gy_ext;
            double gp = r * (-sinp * gx_ext + cosp * gy_ext);

            Fr_pl = M * (gr + r*omega*omega);
            rFp_pl = M * gp;
            Fz_pl = M * gz_ext;

            Fx_pl = M*gx_ext;
            Fy_pl = M*gy_ext;
        }

        /*
         * Just giving the integrals names for ease later.
         */

        double *pl_gas_track = theDomain->pl_gas_track + p*NUM_PL_INTEGRALS;

        //Mass & Momentum
        double dM = pl_gas_track[PL_SNK_M];
        double dPx_grv = pl_gas_track[PL_GRV_PX];
        double dPy_grv = pl_gas_track[PL_GRV_PY];
        double dPz_grv = pl_gas_track[PL_GRV_PZ];
        
        double dPx_snk = pl_gas_track[PL_SNK_PX];
        double dPy_snk = pl_gas_track[PL_SNK_PY];
        double dPz_snk = pl_gas_track[PL_SNK_PZ];

        //Radial & Angular Momentum
        double dJ_grv = pl_gas_track[PL_GRV_JZ];

        double dJ_snk = pl_gas_track[PL_SNK_JZ];
        double dS_snk = pl_gas_track[PL_SNK_SZ];

        //Mass dipole movement
        double dMx_snk = pl_gas_track[PL_SNK_MX];
        double dMy_snk = pl_gas_track[PL_SNK_MY];
        double dMz_snk = pl_gas_track[PL_SNK_MZ];
       
        //Fluid energy added to planet
        double dEf_grv = pl_gas_track[PL_GRV_EGAS];
        double dEf_snk = pl_gas_track[PL_SNK_EGAS];

        //Potential energy btw gas & planet lost during accretion
        double dUpot_snk = pl_gas_track[PL_SNK_UGAS];

        double Uf = pl->Uf;

        /*
         * Calculating some things from the integrals
         */

        //Center of mass movement
        double dx_snk = dMx_snk / M;
        double dy_snk = dMy_snk / M;
        double dz_snk = dMz_snk / M;
        
        double dr_snk = cosp * dx_snk + sinp * dy_snk;
        double dphi_snk = r==0.0 ? 0.0 : (-sinp * dx_snk + cosp * dy_snk) / r;

        // Radial Momentum
        double dPr_grv = cosp * dPx_grv + sinp * dPy_grv;
        double dPr_snk = cosp * dPx_snk + sinp * dPy_snk
                        + M * r*omega * dphi_snk;

        // Orbital Angular Momentum
        double dL_grv = x * dPy_grv - y * dPx_grv;
        double dL_snk = x * dPy_snk - y * dPx_snk 
                        + dMx_snk * vy - dMy_snk * vx;
        
        // Kinetic energy added to the planet
        double dK_grv = vx * dPx_grv + vy * dPy_grv + vz * dPz_grv;
        double dK_snk = vx * dPx_snk + vy * dPy_snk + vz * dPz_snk 
                            - 0.5*v2 * dM;

        // Potential energy added to planet
        double dU_grv = dEf_grv - dK_grv;

        double dUf_snk = - dUpot_snk + (Uf/M) * dM
                        - (dPx_grv/M) * dMx_snk
                        - (dPy_grv/M) * dMy_snk
                        - (dPz_grv/M) * dMz_snk;
        double dUext_snk = Phi_ext * dM
                        - gx_ext * dMx_snk
                        - gy_ext * dMy_snk
                        - gz_ext * dMz_snk;


        double dU_snk = dUf_snk + dUext_snk;
        
        // Internal energy added to planet
        double dEint_snk = dEf_snk - dK_snk - dU_snk;

        double *pl_kin = theDomain->pl_kin + p*NUM_PL_KIN;
        double *pl_aux = theDomain->pl_aux + p*NUM_PL_AUX;


        pl_kin[PL_M] += dM;

        pl_kin[PL_R] += vr_pl * dt + dr_snk;
        pl_kin[PL_PHI] += om_pl * dt + dphi_snk;
        pl_kin[PL_Z] += vz_pl * dt + dz_snk;

        pl_kin[PL_PR] += Fr_pl * dt + dPr_grv + dPr_snk;
        pl_kin[PL_LL] += rFp_pl * dt + dJ_grv + (dJ_snk - dS_snk);
        pl_kin[PL_PZ] += Fz_pl * dt + dPz_grv + dPz_snk;

        pl_kin[PL_SZ] += dS_snk;
        pl_kin[PL_EINT] += dEint_snk;

        //Direct values, easy!
        pl_aux[PL_SNK_M] += dM;
        pl_aux[PL_GRV_PX] += dPx_grv;
        pl_aux[PL_GRV_PY] += dPy_grv;
        pl_aux[PL_GRV_PZ] += dPz_grv;
        pl_aux[PL_GRV_JZ] += dJ_grv;
        pl_aux[PL_SNK_PX] += dPx_snk;
        pl_aux[PL_SNK_PY] += dPy_snk;
        pl_aux[PL_SNK_PZ] += dPz_snk;
        pl_aux[PL_SNK_JZ] += dJ_snk;
        pl_aux[PL_SNK_SZ] += dS_snk;
        pl_aux[PL_SNK_MX] += dMx_snk;
        pl_aux[PL_SNK_MY] += dMy_snk;
        pl_aux[PL_SNK_MZ] += dMz_snk;
        pl_aux[PL_GRV_EGAS] += dEf_grv;
        pl_aux[PL_SNK_EGAS] += dEf_snk;
        pl_aux[PL_SNK_UGAS] += dUf_snk;
        

        // Some more complicated things.
        pl_aux[PL_GRV_LZ] += dL_grv;
        pl_aux[PL_SNK_LZ] += dL_snk;

        pl_aux[PL_GRV_K] += dK_grv;
        pl_aux[PL_SNK_K] += dK_snk;
        
        pl_aux[PL_GRV_U] += dU_grv;
        pl_aux[PL_SNK_U] += dU_snk;

        pl_aux[PL_SNK_EINT] += dEint_snk;

        // The external things
        pl_aux[PL_EXT_PX] += Fx_pl*dt;
        pl_aux[PL_EXT_PY] += Fy_pl*dt;
        pl_aux[PL_EXT_PZ] += Fz_pl*dt;
        pl_aux[PL_EXT_JZ] += rFp_pl*dt;
        pl_aux[PL_EXT_K] += (vx*Fx_pl + vy*Fy_pl + vz*Fz_pl) * dt;
        pl_aux[PL_EXT_U] -= (vx*Fx_pl + vy*Fy_pl + vz*Fz_pl) * dt;

    }
}

