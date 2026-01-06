
#include "disco.h"
#include "geometry.h"
#include "hydro.h"


double minmod( double a , double b , double c ){
   double m = a;
   if( a*b < 0.0 ) m = 0.0;
   if( fabs(b) < fabs(m) ) m = b;
   if( b*c < 0.0 ) m = 0.0;
   if( fabs(c) < fabs(m) ) m = c;
   return(m);
}

void plm_phi( struct domain * theDomain ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double PLM = theDomain->theParList.PLM;
#if ENABLE_CART_INTERP
   int Cartesian_Interp = theDomain->theParList.Cartesian_Interp;
#endif
   int i,j,k,q;


   for( k=0 ; k<Nz ; ++k ){
      for( j=0 ; j<Nr ; ++j ){
         int jk = j+Nr*k;
#if ENABLE_CART_INTERP
         double z = get_centroid(theDomain->z_kph[k],theDomain->z_kph[k-1],2);
         double r = get_centroid(theDomain->r_jph[j],theDomain->r_jph[j-1],1);
         double x[3] = {r, 0.0, z};
         double w = 0.0;
         if(Cartesian_Interp)
            w = getCartInterpWeight(x);
#endif
         for( i=0 ; i<Np[jk] ; ++i ){
            int im = (i == 0) ? Np[jk]-1 : i-1;
            int ip = (i == Np[jk]-1) ? 0 : i+1;
            struct cell * c  = &(theCells[jk][i]);
            struct cell * cL = &(theCells[jk][im]);
            struct cell * cR = &(theCells[jk][ip]);
            double dpL = cL->dphi;
            double dpC = c->dphi;
            double dpR = cR->dphi;
            for( q=0 ; q<NUM_Q ; ++q ){
               double pL = cL->prim[q];
               double pC = c->prim[q];
               double pR = cR->prim[q];
               double sL = pC - pL;
               sL /= .5*( dpC + dpL );
               double sR = pR - pC;
               sR /= .5*( dpR + dpC );
               double sC = pR - pL;
               sC /= .5*( dpL + dpR ) + dpC;
               c->gradp[q] = minmod( PLM*sL , sC , PLM*sR );
            }
#if ENABLE_CART_INTERP
            if(w > 0.0)
            {
               double xL[3] = {r, -0.5*(dpL+dpC), z}; 
               double xC[3] = {r, 0.0, z}; 
               double xR[3] = {r, 0.5*(dpC+dpR), z}; 
               double cartPrimL[NUM_Q];
               double cartPrimC[NUM_Q];
               double cartPrimR[NUM_Q];
               double cartGrad[NUM_Q];
               double grad[NUM_Q];
               geom_rebase_to_cart(cL->prim, xL, cartPrimL);
               geom_rebase_to_cart(c->prim,  xC, cartPrimC);
               geom_rebase_to_cart(cR->prim, xR, cartPrimR);
               for( q=0 ; q<NUM_Q ; ++q ){
                  double pL = cartPrimL[q];
                  double pC = cartPrimC[q];
                  double pR = cartPrimR[q];
                  double sL = pC - pL;
                  double sR = pR - pC;
                  double sC = pR - pL;
                  if(q == URR || q == UPP)
                  {
                      sL /= sin(.5*( dpC + dpL ));
                      sR /= sin(.5*( dpR + dpC ));
                      sC /= 2*sin(0.5*(.5*( dpL + dpR ) + dpC));
                  }
                  else
                  {
                      sL /= .5*( dpC + dpL );
                      sR /= .5*( dpR + dpC );
                      sC /= .5*( dpL + dpR ) + dpC;
                  }
                  cartGrad[q] = minmod( PLM*sL , sC , PLM*sR );
               }
               geom_gradCart_to_grad(cartGrad, c->prim, xC, grad, 0);
               for( q=0 ; q<NUM_Q ; ++q )
                  c->gradp[q] = (1-w) * c->gradp[q] + w * grad[q];
            }
#endif
         }
      }
   }
}

void plm_geom_boundary(struct domain *theDomain, int jmin, int jmax, int kmin,
                        int kmax, int dim, int LR);

void plm_trans( struct domain * theDomain , struct face * theFaces , int Nf , int dim ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double PLM = theDomain->theParList.PLM;
#if ENABLE_CART_INTERP
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int Cartesian_Interp = theDomain->theParList.Cartesian_Interp;
#endif
   
   int i,j,k,q;

   //Clear gradients
   for( j=0 ; j<Nr ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
             if(dim == 1)
                memset( theCells[jk][i].gradr , 0 , NUM_Q*sizeof(double) );
             else
                memset( theCells[jk][i].gradz , 0 , NUM_Q*sizeof(double) );
         }
      }
   }
#if ENABLE_CART_INTERP
   if(Cartesian_Interp && dim == 1)
   {
      for( k=0 ; k<Nz ; ++k ){
         double z = get_centroid(z_kph[k], z_kph[k-1], 2);
         for( j=0 ; j<Nr ; ++j ){
            int jk = j+Nr*k;
            double r = get_centroid(r_jph[j], r_jph[j-1], 1);
            double x[3] = {r, 0.0, z};
            double w = getCartInterpWeight(x);
            for( i=0 ; i<Np[jk] ; ++i ){
                theCells[jk][i].gradr[UPP] = -w*theCells[jk][i].prim[UPP] / r;
             }
          }
       }
   }
#endif

   if(PLM == 0.0)
       return;

   //printf("Doing grads!\n");
   
    //Calc total weight
    for(j = 0; j < Nr; j++)
        for(k = 0; k < Nz; k++)
        {
            int jk = j+Nr*k;
            for(i = 0; i < Np[jk]; i++)
                theCells[jk][i].tempDoub = 0.0;
        }
    int n;
    for(n = 0; n < Nf; n++)
    {
        struct face *f  = &(theFaces[n]);
        (f->L)->tempDoub += f->dA;
        (f->R)->tempDoub += f->dA;
    }
   
    //Add weighted slopes
    for(n=0; n<Nf; ++n)
    {
        struct face * f  = &( theFaces[n] );
        double phi = f->cm[1];
        struct cell * cL = f->L;
        struct cell * cR = f->R;
        double dxL = f->dxL;
        double dxR = f->dxR;
        double phiL = cL->piph - .5*cL->dphi;
        double phiR = cR->piph - .5*cR->dphi;
        double dpL = get_signed_dp(phi,phiL);
        double dpR = get_signed_dp(phi,phiR);
        double dA  = f->dA;

        double dA_fL = dA / cL->tempDoub;
        double dA_fR = dA / cR->tempDoub;

        double *gradL = dim==1 ? cL->gradr : cL->gradz;
        double *gradR = dim==1 ? cR->gradr : cR->gradz;

#if ENABLE_CART_INTERP
        double w = 0.0;
        if(theDomain->theParList.Cartesian_Interp)
            w = getCartInterpWeight(f->cm);

        if(w > 0.0)
        {
            double gradCIL[NUM_Q], gradCIR[NUM_Q];
            geom_cart_interp_grad_trans(cL->prim, cR->prim,
                                        cL->gradp, cR->gradp,
                                        dpL, dpR, f->cm, -dxL, dxR,
                                        gradCIL, gradCIR, dim);
            for(q=0; q<NUM_Q; q++)
            {
                double WL = cL->prim[q] + dpL*cL->gradp[q];
                double WR = cR->prim[q] + dpR*cR->gradp[q];
                double S = (WR-WL)/(dxR+dxL);

                gradL[q] += ((1-w) * S + w * gradCIL[q]) * dA_fL;
                gradR[q] += ((1-w) * S + w * gradCIR[q]) * dA_fR;
            }
        }
        else
        {
#endif
            for( q=0 ; q<NUM_Q ; ++q )
            {
                double WL = cL->prim[q] + dpL*cL->gradp[q];
                double WR = cR->prim[q] + dpR*cR->gradp[q];
                double S = (WR-WL)/(dxR+dxL);
            
                gradL[q] += S*dA_fL;
                gradR[q] += S*dA_fR;
            }
#if ENABLE_CART_INTERP
        }
#endif
    }

    //Slope Limiting
    
    for( n=0 ; n<Nf ; ++n )
    {
        struct face * f  = &( theFaces[n] );
        double phi = f->cm[1];
        struct cell * cL = f->L;
        struct cell * cR = f->R;
        double dxL = f->dxL;
        double dxR = f->dxR;
        double phiL = cL->piph - .5*cL->dphi;
        double phiR = cR->piph - .5*cR->dphi;
        double dpL = get_signed_dp(phi,phiL);
        double dpR = get_signed_dp(phi,phiR);
      
        double *gradL = dim==1 ? cL->gradr : cL->gradz;
        double *gradR = dim==1 ? cR->gradr : cR->gradz;
      
#if ENABLE_CART_INTERP
        double w = 0.0;
        if(Cartesian_Interp)
            w = getCartInterpWeight(f->cm);
        if(w > 0.0)
        {
            double gradCIL[NUM_Q], gradCIR[NUM_Q];
            geom_cart_interp_grad_trans(cL->prim, cR->prim,
                                        cL->gradp, cR->gradp,
                                        dpL, dpR, f->cm, -dxL, dxR,
                                        gradCIL, gradCIR, dim);

            for(q=0; q<NUM_Q; q++)
            {
                double WL = cL->prim[q] + dpL*cL->gradp[q];
                double WR = cR->prim[q] + dpR*cR->gradp[q];
                double S = (WR-WL)/(dxR+dxL);
                double SL = gradL[q];
                double SR = gradR[q];
                double SfL = (1-w) * S + w * gradCIL[q];
                double SfR = (1-w) * S + w * gradCIR[q];

                if(dim == 1 && q == UPP)
                {
                    double SL0 = -w*cL->prim[UPP] / (f->cm[0] - dxL);
                    double SR0 = -w*cR->prim[UPP] / (f->cm[0] + dxL);

                    SL -= SL0;
                    SR -= SR0;
                    
                    if( SfL*SL < 0.0 )
                        gradL[q] = SL0; 
                    else if( fabs(PLM*SfL) < fabs(SL) )
                        gradL[q] = PLM*SfL + SL0;

                    if( SfR*SR < 0.0 )
                        gradR[q] = SR0; 
                    else if( fabs(PLM*SfR) < fabs(SR) )
                        gradR[q] = PLM*SfR + SR0;
                }
                else
                {
                    if( SfL*SL < 0.0 )
                        gradL[q] = 0.0; 
                    else if( fabs(PLM*SfL) < fabs(SL) )
                        gradL[q] = PLM*SfL;
                    

                    if( SfR*SR < 0.0 )
                        gradR[q] = 0.0; 
                    else if( fabs(PLM*SfR) < fabs(SR) )
                        gradR[q] = PLM*SfR;
                }
            }
        }
        else
        {
#endif
            for( q=0 ; q<NUM_Q ; ++q )
            {
                double WL = cL->prim[q] + dpL*cL->gradp[q];
                double WR = cR->prim[q] + dpR*cR->gradp[q];

                double S = (WR-WL)/(dxR+dxL);
                double SL = gradL[q];
                double SR = gradR[q];
                if( S*SL < 0.0 )
                    gradL[q] = 0.0; 
                else if( fabs(PLM*S) < fabs(SL) )
                    gradL[q] = PLM*S;

                if( S*SR < 0.0 )
                    gradR[q] = 0.0; 
                else if( fabs(PLM*S) < fabs(SR) )
                    gradR[q] = PLM*S;
            }
#if ENABLE_CART_INTERP
        }
#endif
    }
   
   // Geometric boundaries don't have ghost zones, so the gradients there need
   // to be fixed.  
 

   if(dim == 1 && theDomain->NgRa == 0)
      plm_geom_boundary(theDomain, 0, 0, 0, Nz-1, dim, 0);

   if(dim == 1 && theDomain->NgRb == 0)
      plm_geom_boundary(theDomain, Nr-1, Nr-1, 0, Nz-1, dim, 1);

   if(dim == 2 && theDomain->NgZa == 0)
      plm_geom_boundary(theDomain, 0, Nr-1, 0, 0, dim, 0);

   if(dim == 2 && theDomain->NgZb == 0)
      plm_geom_boundary(theDomain, 0, Nr-1, Nz-1, Nz-1, dim, 1);
}

void plm_geom_boundary(struct domain *theDomain, int jmin, int jmax, 
                        int kmin, int kmax, int dim, int LR)
{
    struct cell ** theCells = theDomain->theCells;
    int Nr = theDomain->Nr;
    int * Np = theDomain->Np;
    double * r_jph = theDomain->r_jph;
    double * z_kph = theDomain->z_kph;
    double PLM = theDomain->theParList.PLM;

    int i,j,k;

    double xp[3], xm[3];
    for(k=kmin; k<=kmax; k++)
    {
       xp[2] = z_kph[k];
       xm[2] = z_kph[k-1];
       for(j=jmin; j<=jmax; j++)
       {
           int jk = j+Nr*k;
           xp[0] = r_jph[j];
           xm[0] = r_jph[j-1];
           for(i=0; i < Np[jk]; i++)
           {
               struct cell * c = &(theCells[jk][i]);
               xp[1] = c->piph;
               xm[1] = c->piph - c->dphi;

               if(dim == 1)
                   geom_grad(c->prim, c->gradr, xp, xm, PLM, dim, LR); 
               else
                   geom_grad(c->prim, c->gradz, xp, xm, PLM, dim, LR); 
           }
       }
    }
}
