enum{_HLL_,_HLLC_,_HLLD_,_HLLC_DAMPCENTER_};

#include "paul.h"
#include "hydro.h"
#include "geometry.h"
#include "riemann.h"

static int mesh_motion = 0;
static int riemann_solver = 0;
static int visc_flag = 0;
static int use_B_fields = 0;
#if ENABLE_CART_INTERP
static int Cartesian_Interp = 0;
#endif


void setRiemannParams( struct domain * theDomain ){
   mesh_motion = theDomain->theParList.Mesh_Motion;
   riemann_solver = theDomain->theParList.Riemann_Solver;
   visc_flag = theDomain->theParList.visc_flag;
   use_B_fields = set_B_flag();
#if ENABLE_CART_INTERP
   Cartesian_Interp = theDomain->theParList.Cartesian_Interp;
#endif

   if( !use_B_fields && riemann_solver == _HLLD_ && theDomain->rank==0 ){
      printf("Ya dun goofed.\nRiemann Solver = HLLD,\nHydro does not include magnetic fields.\n");
      exit(1);
   }
   if( !use_B_fields && (NUM_C > 5 || NUM_EDGES > 0 || NUM_FACES > 0 || NUM_AZ_EDGES > 0 ) && theDomain->rank==0 ){
      printf("Warning:  You are not solving MHD equations but possibly storing more variables than you need.\nNum Conserved Vars = %d\nFaces = %d\nEdges = %d\nAzimuthal Edges = %d\nCode might still work fine, this is just a warning.\n",NUM_C,NUM_FACES,NUM_EDGES,NUM_AZ_EDGES);
   }

}

void get_Ustar_HLLD( double , const double * , const double * , 
                    double * , double * , const double * , const double * );

void riemann_phi( struct cell * cL , struct cell * cR, double * x ,
                const double *xp, const double *xm, double dAdt ){

   double primL[NUM_Q];
   double primR[NUM_Q];
   double gradp[NUM_Q];

   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      primL[q] = cL->prim[q] + .5*cL->gradp[q]*cL->dphi;
      primR[q] = cR->prim[q] - .5*cR->gradp[q]*cR->dphi;
      gradp[q] = (cR->prim[q] - cL->prim[q]) / (0.5*(cR->dphi + cL->dphi));
   }

   double xc[3] = {x[0], x[1], x[2]};

#if ENABLE_CART_INTERP
   if(Cartesian_Interp)
   {
       double weight = getCartInterpWeight(x);
       
       if(weight > 0.0)
       {
           double xL[3] = {x[0], cL->piph - 0.5*cL->dphi, x[2]};
           double xR[3] = {x[0], cR->piph - 0.5*cR->dphi, x[2]};
           xc[0] = 0.5*(xp[0] + xm[0]);
           double dr = xc[0] - x[0];
           geom_interpolate(cL->prim, cL->gradr, cL->gradp, cL->gradz,
                            xL,  dr, 0.5*cL->dphi, 0.0, primL, weight);
           geom_interpolate(cR->prim, cR->gradr, cR->gradp, cR->gradz,
                            xR,  dr, -0.5*cR->dphi, 0.0, primR, weight);
       }
   }
#endif


   double n[3] = {0.0,1.0,0.0};
   double hn = get_scale_factor(xc, 0);

   if( use_B_fields && NUM_Q > BPP ){
      double Bp = .5*(primL[BPP]+primR[BPP]);
      primL[BPP] = Bp;
      primR[BPP] = Bp;
   }

   double Er,Ez,Br,Bz;
   solve_riemann(primL , primR , cL->cons , cR->cons ,
                 cL->gradr, gradp, cL->gradz,
                 cR->gradr, gradp, cR->gradz,
                 xc , n , xp, xm, hn*cL->wiph , dAdt , 0 ,
                 &Ez , &Br , &Er , &Bz , NULL, NULL);
   /*
   solve_riemann(primL , primR , cL->cons , cR->cons ,
                 cL->gradr, cL->gradp, cL->gradz,
                 cR->gradr, cR->gradp, cR->gradz,
                 x , n , hn*cL->wiph , dAdt , 0 , &Ez , &Br , &Er , &Bz );
    */

#if NUM_EDGES == 4
    cL->E[0] = .5*Ez;
    cL->E[1] = .5*Ez;

    cR->E[2] = .5*Ez;
    cR->E[3] = .5*Ez;

    cL->B[0] = .5*Br;
    cL->B[1] = .5*Br;

    cR->B[2] = .5*Br;
    cR->B[3] = .5*Br;
#elif NUM_EDGES == 8
    cL->E[0] = .5*Ez;
    cL->E[1] = .5*Ez;

    cR->E[2] = .5*Ez;
    cR->E[3] = .5*Ez;

    cL->B[0] = .5*Br;
    cL->B[1] = .5*Br;

    cR->B[2] = .5*Br;
    cR->B[3] = .5*Br;

    cL->E[4] = .5*Er;
    cL->E[5] = .5*Er;

    cR->E[6] = .5*Er;
    cR->E[7] = .5*Er;

    cL->B[4] = .5*Bz;
    cL->B[5] = .5*Bz;

    cR->B[6] = .5*Bz;
    cR->B[7] = .5*Bz;
#endif
}

void riemann_trans( struct face * F , double dt , int dim , double rp,
                    double rm, double zp, double zm,
                    double *fdAdt_hydro, double *fdAdt_visc){

   struct cell * cL = F->L;
   struct cell * cR = F->R;
   double dAdt      = F->dA*dt;
   double dxL       = F->dxL;
   double dxR       = F->dxR;
   double phi       = F->cm[1];
   double dphi      = F->dphi;

   double xp[3] = {rp, phi+0.5*dphi, zp};
   double xm[3] = {rm, phi-0.5*dphi, zm};

   double primL[NUM_Q];
   double primR[NUM_Q];
   double grad[NUM_Q];

   double phiL = cL->piph - .5*cL->dphi;
   double phiR = cR->piph - .5*cR->dphi;
   double dpL = get_signed_dp(phi,phiL);
   double dpR = get_signed_dp(phi,phiR);
      
   double *gradL = dim==1 ? cL->gradr : cL->gradz;
   double *gradR = dim==1 ? cR->gradr : cR->gradz;

   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      primL[q] = cL->prim[q] + gradL[q]*dxL + cL->gradp[q]*dpL;
      primR[q] = cR->prim[q] - gradR[q]*dxR + cR->gradp[q]*dpR;
      grad[q] = ((cR->prim[q] + cR->gradp[q]*dpR)
                  - (cL->prim[q] + cL->gradp[q]*dpL)) / (dxL+dxR);
   }
    
   double x[3] = {F->cm[0], F->cm[1], F->cm[2]};

#if ENABLE_CART_INTERP
   if(Cartesian_Interp)
   {
       double weight = getCartInterpWeight(F->cm);

       if(weight > 0.0)
       {
           double xL[3] = {F->cm[0], phiL, F->cm[2]};
           double xR[3] = {F->cm[0], phiR, F->cm[2]};
           double drL = 0.0;
           double drR = 0.0;
           double dzL = 0.0;
           double dzR = 0.0;
           if(dim == 1)
           {
               xL[0] -= dxL;
               xR[0] += dxR;
               drL = dxL;
               drR = -dxR;
               x[0] *= cos(0.5*F->dphi);
               drL = x[0] - xL[0];
               drR = x[0] - xR[0];
           }
           else
           {
               xL[2] -= dxL;
               xR[2] += dxR;
               dzL = dxL;
               dzR = -dxR;
           }
           geom_interpolate(cL->prim, cL->gradr, cL->gradp, cL->gradz,
                            xL, drL, dpL, dzL, primL, weight);
           geom_interpolate(cR->prim, cR->gradr, cR->gradp, cR->gradz,
                            xR, drR, dpR, dzR, primR, weight);
       }
   }
#endif

   double n[3] = {0.0,0.0,0.0};
   if( dim==1 ) n[0] = 1.0; else n[2] = 1.0;

   if( use_B_fields && NUM_C > BZZ){
      int BTRANS;
      if( dim==1 ) BTRANS = BRR; else BTRANS = BZZ;
      double Bavg = .5*(primL[BTRANS]+primR[BTRANS]);
      primL[BTRANS] = Bavg;
      primR[BTRANS] = Bavg;
   }

   double Erz,Brz,Ephi;
    /*
   solve_riemann(primL , primR , cL->cons , cR->cons ,
                 cL->gradr, cL->gradp, cL->gradz,
                 cR->gradr, cR->gradp, cR->gradz,
                 F->cm , n , 0.0 , dAdt , dim , &Erz , &Brz , &Ephi , NULL );
    */
   if(dim == 1)
      solve_riemann(primL , primR , cL->cons , cR->cons ,
                    grad, cL->gradp, cL->gradz,
                    grad, cR->gradp, cR->gradz,
                    x , n , xp, xm, 0.0 , dAdt , dim ,
                    &Erz , &Brz , &Ephi, NULL ,
                    fdAdt_hydro, fdAdt_visc);
   else
      solve_riemann(primL , primR , cL->cons , cR->cons ,
                    cL->gradr, cL->gradp, grad,
                    cR->gradr, cR->gradp, grad,
                    x , n , xp, xm, 0.0 , dAdt , dim ,
                    &Erz , &Brz , &Ephi, NULL,
                    fdAdt_hydro, fdAdt_visc);
 

#if NUM_EDGES >= 4
   double fracL = F->dphi / cL->dphi;
   double fracR = F->dphi / cR->dphi;
   if(dim == 1) {
      cL->E[1] += .5*Erz*fracL;
      cL->E[3] += .5*Erz*fracL;

      cR->E[0] += .5*Erz*fracR;
      cR->E[2] += .5*Erz*fracR;

      cL->B[1] += .5*Brz*fracL;
      cL->B[3] += .5*Brz*fracL;

      cR->B[0] += .5*Brz*fracR;
      cR->B[2] += .5*Brz*fracR;
   }
#if NUM_EDGES == 8
   if(dim==2){
      cL->E[5] += .5*Erz*fracL;
      cL->E[7] += .5*Erz*fracL;

      cR->E[4] += .5*Erz*fracR;
      cR->E[6] += .5*Erz*fracR;

      cL->B[5] += .5*Brz*fracL;
      cL->B[7] += .5*Brz*fracL;

      cR->B[4] += .5*Brz*fracR;
      cR->B[6] += .5*Brz*fracR;
   }
#endif
#endif
#if NUM_AZ_EDGES == 4
   if(dim==1 ){
      if(F->LRtype==0)
         cL->E_phi[1] = Ephi;
      else
         cR->E_phi[0] = Ephi;
   }
   else if( dim==2 ){
      if(F->LRtype==0)
         cL->E_phi[3] = Ephi;
      else
         cR->E_phi[2] = Ephi;
   }
#endif
}


void solve_riemann(const double *primL, const double *primR,
                   double *consL, double *consR,
                   const double *gradLr, const double *gradLp,
                   const double *gradLz,
                   const double *gradRr, const double *gradRp, 
                   const double *gradRz,
                   const double *x, const double *n, 
                   const double *xp, const double *xm,
                   double w, double dAdt, 
                   int dim,
                   double *E1_riemann, double *B1_riemann,
                   double *E2_riemann, double *B2_riemann,
                   double *fdAdt_hydro, double *fdAdt_visc)
{

   int q;

   double Flux[NUM_Q];
   double Ustr[NUM_Q];
   for( q=0 ; q<NUM_Q ; ++q ){
      Flux[q] = 0.0;
      Ustr[q] = 0.0;
   }

   double r = x[0];

   if( riemann_solver == _HLL_ || riemann_solver == _HLLC_ || riemann_solver == _HLLC_DAMPCENTER_ ){
      double Sl,Sr,Ss;
      double Bpack[5];
      vel( primL , primR , &Sl , &Sr , &Ss , n , x , Bpack );

      if( w < Sl ){
         flux( primL , Flux , x , n , xp , xm);
         prim2cons( primL , Ustr , x , 1.0, NULL, NULL);

      }else if( w > Sr ){
         flux( primR , Flux , x , n , xp , xm);
         prim2cons( primR , Ustr , x , 1.0, NULL, NULL);

      }else{
         if( riemann_solver == _HLL_ || (r < 0.1 && riemann_solver == _HLLC_DAMPCENTER_) ){
            double Fl[NUM_Q];
            double Fr[NUM_Q];
            double Ul[NUM_Q];
            double Ur[NUM_Q];
   
            double aL =  Sr;
            double aR = -Sl;
 
            prim2cons( primL , Ul , x , 1.0 , NULL, NULL);
            prim2cons( primR , Ur , x , 1.0 , NULL, NULL);
            flux( primL , Fl , x , n , xp , xm);
            flux( primR , Fr , x , n , xp , xm);

            for( q=0 ; q<NUM_Q ; ++q ){
               Flux[q] = ( aL*Fl[q] + aR*Fr[q] + aL*aR*( Ul[q] - Ur[q] ) )/( aL + aR );
               Ustr[q] = ( aR*Ul[q] + aL*Ur[q] + Fl[q] - Fr[q] )/( aL + aR );
            }
         }else{
            double Uk[NUM_Q];
            double Fk[NUM_Q];
            if( w < Ss ){
               prim2cons( primL , Uk , x , 1.0, NULL, NULL);
               getUstar( primL , Ustr , x , Sl , Ss , n , Bpack ); 
               flux( primL , Fk , x , n , xp , xm); 

               for( q=0 ; q<NUM_Q ; ++q ){
                  Flux[q] = Fk[q] + Sl*( Ustr[q] - Uk[q] );
               }    
            }else{
               prim2cons( primR , Uk , x , 1.0, NULL, NULL);
               getUstar( primR , Ustr , x , Sr , Ss , n , Bpack ); 
               flux( primR , Fk , x , n , xp , xm); 

               for( q=0 ; q<NUM_Q ; ++q ){
                  Flux[q] = Fk[q] + Sr*( Ustr[q] - Uk[q] );
               } 
            } 
         }
      }
   }else{
      get_Ustar_HLLD( w , primL , primR , Flux , Ustr , x , n );
   }

   double vFlux[NUM_Q] = {0};
   if( visc_flag ){
      double prim[NUM_Q];
      double gradr[NUM_Q];
      double gradp[NUM_Q];
      double gradz[NUM_Q];
      for( q=0 ; q<NUM_Q ; ++q ){
         prim[q] = .5*(primL[q]+primR[q]);
         gradr[q] = .5*(gradLr[q]+gradRr[q]);
         gradp[q] = .5*(gradLp[q]+gradRp[q]);
         gradz[q] = .5*(gradLz[q]+gradRz[q]);
         vFlux[q] = 0.0;
      }
      visc_flux(prim, gradr, gradp, gradz, vFlux, x, n);
   }

   /*
   FILE *f = fopen("flux.dat", "a");
    fprintf(f, "%d %.15e %.15e %.15e", dim, x[0], x[1], x[2]);
   fprintf(f, " %.15e %.15e %.15e", xp[0], xp[1], xp[2]);
   fprintf(f, " %.15e %.15e %.15e", xm[0], xm[1], xm[2]);
   fprintf(f, " %.15e", dAdt);
   for(q=0; q<NUM_Q; q++) fprintf(f, " %.15e", primL[q]);
   for(q=0; q<NUM_Q; q++) fprintf(f, " %.15e", primR[q]);
   for(q=0; q<NUM_Q; q++) fprintf(f, " %.15e", Flux[q]);
   fprintf(f, "\n");
   fclose(f);
   */

   for( q=0 ; q<NUM_Q ; ++q ){
      consL[q] -= (Flux[q] + vFlux[q] - w*Ustr[q])*dAdt;
      consR[q] += (Flux[q] + vFlux[q] - w*Ustr[q])*dAdt;
   }

   if(fdAdt_hydro != NULL)
   {
      for( q=0 ; q<NUM_Q ; ++q )
         fdAdt_hydro[q] = (Flux[q] - w*Ustr[q])*dAdt;
   }

   if(fdAdt_visc != NULL)
   {
      for( q=0 ; q<NUM_Q ; ++q )
         fdAdt_visc[q] = vFlux[q]*dAdt;
   }

   flux_to_E( Flux , Ustr , x , E1_riemann , B1_riemann , E2_riemann , B2_riemann , dim );
}

