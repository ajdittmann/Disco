#include <stdio.h>
#include <math.h>
#include "paul.h"
#include "Hydro/metric.h"
#include "Hydro/frame.h"

#define ACC 1e-12
#define MAX_ITER 10

#define DEBUG 0
#define DEBUG2 0
#define DEBUG3 0
#define DEBUG4 0

enum{DD,EN,MX,MY,MZ,BX,BY,BZ};

static double gamma_law = 0.0;
static int isothermal = 0;
static double RMIN = 0.0;
static double RMAX = 0.0;

//From GRMHD
void prim2cons( double * , double * , double * , double );
void flux( double * , double * , double * , double * );
void vel( double * , double * , double * , double * , double * , double * , double * , double * );
int calc_msfast(double *velL, double *velR, double rhoh, double cs2, double w, 
                double vi, double b2, double b0, double bi, double al, 
                double bei, double igamii, double *x);

//Global Functions
void get_Ustar_HLLD(double w, double *pL, double *pR, double *F, double *U, 
                    double *x, double *n);
double get_cs2( double );

//Local Functions
void calc_tetrads(double al, double *be, double *gam, double *igam, 
                    double sqrtgam, int xx, double *e, double *ie);
void calc_UF(double rho, double P, double *u, double *B, double *U, double *F);
int solve_HLLD_SR(double sL, double sR, double Bx, double *UL, double *FL,
                    double *UR, double *FR, double w, double *pL, double *pR,
                    double *U, double *F);
int checkSolution(double Ps, double sL, double sR, double *vaL, double *vaR, 
                    double *BaL, double *BaR, double enthL, double enthR, 
                    double *KaL, double *KaR, double *Bc, double *vcL, 
                    double *vcR);
void calc_va(double p, double s, double Bx, double *R, double *va, 
                double *dva);
void calc_Ba(double s, double Bx, double *R, double *va, double *Ba, 
                double *dva, double *dBa);
void calc_entha(double p, double s, double *R, double *va, double *enth,
                    double *dva, double *denth);
void calc_Ka(double p, double s, double Bx, double enth, double *R, double *Ka,
                double denth, double *dKa, int LR);
void cont_func(double Bx, double *BaL, double *BaR,
                double *vaL, double *vaR, double *KaL, double *KaR,
                double enthL, double enthR, double *f,
                double *dBaL, double *dBaR, double *dvaL, double *dvaR,
                double *dKaL, double *dKaR, double denthL, double denthR,
                double *df);
void calc_Bc(double Bx, double *BaL, double *BaR, double *vaL, double *vaR,
                double saL, double saR, double *Bc);
void calc_vc(double Bx, double *Ka, double *Bc, double enth, double *vc,
                int LR);
int calc_Ps(double *RL, double *RR, double Bx, double sL, double sR, 
                double *Ps);
void calc_state(double p, double *RL, double *RR, double Bx, double sL, 
                    double sR, double *vaL, double *vaR, double *BaL, 
                    double *BaR, double *enthL, double *enthR, double *KaL,
                    double *KaR, double *Bc, double *vcL, double *vcR);
void calc_Ua(double p, double s, double Bx, double *R, double *va, double *Ba,
                double *Ua);
void calc_Uc(double p, double sa, double Bx, double *va, double *vc, 
                double *Bc, double *Ua, double *Uc);
double cons2prim_hlld(double *U, double *p);
double cons2prim_noble2d_hlld(double *U, double *p);
double cons2prim_isothermal_hlld(double *U, double *p);

void setHlldParams( struct domain * theDomain )
{
   gamma_law = theDomain->theParList.Adiabatic_Index;
   isothermal = theDomain->theParList.isothermal_flag;
   double rmin = theDomain->theParList.rmin;
   double rmax = theDomain->theParList.rmax;
   int numr = theDomain->theParList.Num_R;
   double dr = (rmax-rmin) / numr;
   RMIN = rmin + 10*dr;
   RMAX = rmax - 10*dr;
}

void get_Ustar_HLLD(double w, double *pL, double *pR, double *F, double *U, 
                    double *x, double *n)
{
    if(DEBUG2)
        printf("In get_Ustar_HLLD\n");
    int i,j,q;
    double r = x[0];
    double al, be[3], gam[9], igam[9], jac;
    al = metric_lapse(x);
    metric_shift(x, be);
    metric_gam(x, gam);
    metric_igam(x, igam);
    jac = metric_jacobian(x) / r;
    double sqrtgam = jac / al;

    double rhoL = pL[RHO];
    double PL = pL[PPP];
    double lL[3] = {pL[URR], pL[UPP], pL[UZZ]};
    double BL[3] = {pL[BRR]/sqrtgam,pL[BPP]/(r*sqrtgam),pL[BZZ]/sqrtgam};

    double rhoR = pR[RHO];
    double PR = pR[PPP];
    double lR[3] = {pR[URR], pR[UPP], pR[UZZ]};
    double BR[3] = {pR[BRR]/sqrtgam,pR[BPP]/(r*sqrtgam),pR[BZZ]/sqrtgam};

    double Uf[4];
    frame_U(x, Uf);

    // Figure out direction.
    int xx;
    if(n[0] > n[1] && n[0] > n[2])
        xx = 0;
    else if(n[1] > n[2])
        xx = 1;
    else
        xx = 2;
    
    //Face 4-velocity
    if(xx == 1)
        w /= r;
    double uw[4];
    double be2 = 0.0;
    double bex = 0.0;
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            be2 += gam[3*i+j]*be[i]*be[j];
    for(i=0; i<3; i++)
        bex += gam[3*xx+i]*be[i];
    uw[0] = 1.0/sqrt(al*al - be2 - 2*bex*w - gam[3*xx+xx]*w*w);
    uw[1] = 0.0;
    uw[2] = 0.0;
    uw[3] = 0.0;
    uw[xx+1] = uw[0]*w;

    if(DEBUG4)
    {
        printf("al = %.8lg, be = (%.8lg %.8lg %.8lg)\n", al, 
                    be[0], be[1], be[2]);
        printf("be2 = %.8lg, bex = %.8lg, gamxx = %.8lg\n", 
                    be2, bex, gam[3*xx+xx]);
        printf("w = %.8lg uw = (%.8lg, %.8lg, %.8lg, %.8lg)\n", 
                w, uw[0], uw[1], uw[2], uw[3]);
    }

    //Calculate orthonormal transform matrices
    double e[16], ie[16];
    calc_tetrads(al, be, gam, igam, r*sqrtgam, xx, e, ie);

    //Transform vector prims to orthonormal frame
    //RELIES on e_(t)^\mu = n^\mu
    
    double uLo[4], uRo[4], BLo[3], BRo[3], uwo[4];

    for(i=0; i<3; i++)
    {
        uLo[i+1] = 0.0;
        uRo[i+1] = 0.0;
        BLo[i] = 0.0;
        BRo[i] = 0.0;
        for(j=0; j<3; j++)
        {
            uLo[i+1] += e[4*(i+1)+j+1]*lL[j];
            uRo[i+1] += e[4*(i+1)+j+1]*lR[j];
            BLo[i] += ie[4*(i+1)+j+1]*BL[j];
            BRo[i] += ie[4*(i+1)+j+1]*BR[j];
        }
    }
    uLo[0] = sqrt(1.0 + uLo[1]*uLo[1] + uLo[2]*uLo[2] + uLo[3]*uLo[3]);
    uRo[0] = sqrt(1.0 + uRo[1]*uRo[1] + uRo[2]*uRo[2] + uRo[3]*uRo[3]);

    uwo[0] = ie[4*0+0]*uw[0] + ie[4*0+xx+1]*uw[xx+1];
    uwo[1] = ie[4*1+0]*uw[0] + ie[4*1+xx+1]*uw[xx+1];
    uwo[2] = ie[4*2+0]*uw[0] + ie[4*2+xx+1]*uw[xx+1];
    uwo[3] = ie[4*3+0]*uw[0] + ie[4*3+xx+1]*uw[xx+1];

    //TODO: check this is kosher? Bx should already satisfy this...
    double Bx = 0.5*(BLo[0] + BRo[0]);
    BLo[0] = Bx;
    BRo[0] = Bx;
    
    //Transform wave/face speeds sL,sR,w to orthonormal frame.
    //TODO: Verify correctness.
    double wo = uwo[1]/uwo[0];
    double sLo, sRo, sL1, sL2, sR1, sR2;

    double rhohL = rhoL + gamma_law/(gamma_law-1)*PL;
    double rhohR = rhoR + gamma_law/(gamma_law-1)*PR;
    double cs2L = gamma_law * PL / rhohL;
    double cs2R = gamma_law * PR / rhohR;
    double b0Lo = uLo[1]*BLo[0] + uLo[2]*BLo[1] + uLo[3]*BLo[2];
    double b0Ro = uRo[1]*BRo[0] + uRo[2]*BRo[1] + uRo[3]*BRo[2];
    double b2L = (BLo[0]*BLo[0]+BLo[1]*BLo[1]+BLo[2]*BLo[2] + b0Lo*b0Lo)
                    / (uLo[0]*uLo[0]);
    double b2R = (BRo[0]*BRo[0]+BRo[1]*BRo[1]+BRo[2]*BRo[2] + b0Ro*b0Ro)
                    / (uRo[0]*uRo[0]);
    double bxLo = (BLo[0] + b0Lo*uLo[1]) / uLo[0];
    double bxRo = (BLo[0] + b0Ro*uRo[1]) / uRo[0];

    int err = calc_msfast(&sL1, &sL2, rhohL, cs2L, uLo[0], uLo[1]/uLo[0], 
                            b2L, b0Lo, bxLo, 1.0, 0.0, 1.0, x);
    err = calc_msfast(&sR1, &sR2, rhohR, cs2R, uRo[0], uRo[1]/uRo[0], 
                            b2R, b0Ro, bxRo, 1.0, 0.0, 1.0, x);
    sLo = sL1 < sR1 ? sL1 : sR1;
    sRo = sL2 > sR2 ? sL2 : sR2;

    //Calculate U,F L/R in orthonormal frame.
    double ULo[NUM_C], FLo[NUM_C], URo[NUM_C], FRo[NUM_C], 
            Uo[NUM_C], Fo[NUM_C];
    calc_UF(rhoL, PL, uLo, BLo, ULo, FLo);
    calc_UF(rhoR, PR, uRo, BRo, URo, FRo);

    //Solve Riemann Problem!
    double primLo[5], primRo[5];
    primLo[RHO] = rhoL;
    primLo[PPP] = PL;
    primLo[URR] = uLo[1];
    primLo[UPP] = uLo[2];
    primLo[UZZ] = uLo[3];
    primRo[RHO] = rhoR;
    primRo[PPP] = PR;
    primRo[URR] = uRo[1];
    primRo[UPP] = uRo[2];
    primRo[UZZ] = uRo[3];
    int zone = solve_HLLD_SR(sLo, sRo, Bx, ULo, FLo, URo, FRo, wo, 
                                primLo, primRo, Uo, Fo);

    if(zone < 0 && x[0]>RMIN && x[0]<RMAX)
    {
        printf("HLLD Failure. Using HLL.  x=(%.6lg %.6lg %.6lg) n=(%.0lf %.0lf %.0lf)\n", 
                x[0], x[1], x[2], n[0], n[1], n[2]);
        printf("%.3lg %.3lg\n", RMIN, RMAX);
    }

    //Convert momenta and energy to coordinate frame
    //TODO: Check accuracy

    double UoTAU, UoSRR, UoLLL, UoSZZ, FoTAU, FoSRR, FoLLL, FoSZZ;
    double UoBRR, UoBPP, UoBZZ, FoBRR, FoBPP, FoBZZ;
    UoTAU = -ie[4*0+0]*Uo[EN] + ie[4*1+0]*Uo[MX]+ ie[4*2+0]*Uo[MY]
           + ie[4*3+0]*Uo[MZ];
    UoSRR = -ie[4*0+1]*Uo[EN] + ie[4*1+1]*Uo[MX]+ ie[4*2+1]*Uo[MY]
           + ie[4*3+1]*Uo[MZ];
    UoLLL = -ie[4*0+2]*Uo[EN] + ie[4*1+2]*Uo[MX]+ ie[4*2+2]*Uo[MY]
           + ie[4*3+2]*Uo[MZ];
    UoSZZ = -ie[4*0+3]*Uo[EN] + ie[4*1+3]*Uo[MX]+ ie[4*2+3]*Uo[MY]
           + ie[4*3+3]*Uo[MZ];
    UoBRR = e[4*1+1]*Uo[BX] + e[4*2+1]*Uo[BY] + e[4*3+1]*Uo[BZ];
    UoBPP = e[4*1+2]*Uo[BX] + e[4*2+2]*Uo[BY] + e[4*3+2]*Uo[BZ];
    UoBZZ = e[4*1+3]*Uo[BX] + e[4*2+3]*Uo[BY] + e[4*3+3]*Uo[BZ];
    FoTAU = -ie[4*0+0]*Fo[EN] + ie[4*1+0]*Fo[MX]+ ie[4*2+0]*Fo[MY]
           + ie[4*3+0]*Fo[MZ];
    FoSRR = -ie[4*0+1]*Fo[EN] + ie[4*1+1]*Fo[MX]+ ie[4*2+1]*Fo[MY]
           + ie[4*3+1]*Fo[MZ];
    FoLLL = -ie[4*0+2]*Fo[EN] + ie[4*1+2]*Fo[MX]+ ie[4*2+2]*Fo[MY]
           + ie[4*3+2]*Fo[MZ];
    FoSZZ = -ie[4*0+3]*Fo[EN] + ie[4*1+3]*Fo[MX]+ ie[4*2+3]*Fo[MY]
           + ie[4*3+3]*Fo[MZ];
    FoBRR = e[4*1+1]*Fo[BX] + e[4*2+1]*Fo[BY] + e[4*3+1]*Fo[BZ];
    FoBPP = e[4*1+2]*Fo[BX] + e[4*2+2]*Fo[BY] + e[4*3+2]*Fo[BZ];
    FoBZZ = e[4*1+3]*Fo[BX] + e[4*2+3]*Fo[BY] + e[4*3+3]*Fo[BZ];

    //Actual fluxes now.
    double hn = 1.0*n[0] + r*n[1] + 1.0*n[2];

    U[DDD] = jac *      (e[4*0+0]*Uo[DD] + e[4*1+0]*Fo[DD]);
    U[SRR] = jac *      (e[4*0+0]*UoSRR  + e[4*1+0]*FoSRR);
    U[LLL] = jac *      (e[4*0+0]*UoLLL  + e[4*1+0]*FoLLL);
    U[SZZ] = jac *      (e[4*0+0]*UoSZZ  + e[4*1+0]*FoSZZ);
    U[TAU] = jac *      (e[4*0+0]*UoTAU  + e[4*1+0]*FoTAU);
    U[BRR] = sqrtgam *      (e[4*0+0]*UoBRR  + e[4*1+0]*FoBRR);
    U[BPP] = sqrtgam *      (e[4*0+0]*UoBPP  + e[4*1+0]*FoBPP);
    U[BZZ] = sqrtgam *      (e[4*0+0]*UoBZZ  + e[4*1+0]*FoBZZ);
    F[DDD] = jac * hn * (e[4*0+xx+1]*Uo[DD] + e[4*1+xx+1]*Fo[DD]);
    F[SRR] = jac * hn * (e[4*0+xx+1]*UoSRR  + e[4*1+xx+1]*FoSRR);
    F[LLL] = jac * hn * (e[4*0+xx+1]*UoLLL  + e[4*1+xx+1]*FoLLL);
    F[SZZ] = jac * hn * (e[4*0+xx+1]*UoSZZ  + e[4*1+xx+1]*FoSZZ);
    F[TAU] = jac * hn * (e[4*0+xx+1]*UoTAU  + e[4*1+xx+1]*FoTAU);
    F[BRR] = sqrtgam * hn * (e[4*0+xx+1]*UoBRR  + e[4*1+xx+1]*FoBRR);
    F[BPP] = sqrtgam * hn * (e[4*0+xx+1]*UoBPP  + e[4*1+xx+1]*FoBPP);
    F[BZZ] = sqrtgam * hn * (e[4*0+xx+1]*UoBZZ  + e[4*1+xx+1]*FoBZZ);

    //Subtract frame velocity and rest-mass density.
    U[TAU] = -Uf[0]*U[TAU]-Uf[1]*U[SRR]-Uf[2]*U[LLL]-Uf[3]*U[SZZ] - U[DDD];
    F[TAU] = -Uf[0]*F[TAU]-Uf[1]*F[SRR]-Uf[2]*F[LLL]-Uf[3]*F[SZZ] - F[DDD];

    //Magnetic Field adjustment
    U[BRR] /= r;
    F[BRR] /= r;

    //Passive scalar.
    for(q=NUM_C; q<NUM_Q; q++)
    {
        if(zone == 0 || zone == 1 || zone == 2)
        {
            U[q] = pL[q]*U[DDD];
            F[q] = pL[q]*F[DDD];
        }
        else
        {
            U[q] = pR[q]*U[DDD];
            F[q] = pR[q]*F[DDD];
        }
    }
    if(DEBUG4)
    {
        double UL[NUM_Q], UR[NUM_Q], FL[NUM_Q], FR[NUM_Q], 
                Uhll[NUM_C], Fhll[NUM_C];
        double sL, sR, sS, Bpack[5];
        prim2cons(pL, UL, x, 1.0);
        prim2cons(pR, UR, x, 1.0);
        flux(pL, FL, x, n);
        flux(pR, FR, x, n);
        vel(pL, pR, &sL, &sR, &sS, n, x, Bpack);
        for(q = 0; q<NUM_C; q++)
        {
            Uhll[q] = (sR*UR[q] - sL*UL[q] + FL[q]-FR[q]) / (sR - sL);
            Fhll[q] = (sR*FL[q] - sL*FR[q] + sL*sR*(UR[q]-UL[q])) / (sR - sL);
        }

        printf("FLUX. n=(%.1lf %.1lf %.1lf) xx=%d sL=%.12lg sR=%.12lg\n",
                    n[0], n[1], n[2], xx, sL, sR);
        printf("      x=(%.12lg %.12lg %.12lg)\n", x[0], x[1], x[2]);
        printf("      primL=(%.12lg %.12lg %.12lg %.12lg)\n", 
                                pL[0], pL[1], pL[2], pL[3]);
        printf("            (%.12lg %.12lg %.12lg %.12lg)\n", 
                                pL[4], pL[5], pL[6], pL[7]);
        printf("      primR=(%.12lg %.12lg %.12lg %.12lg)\n", 
                                pR[0], pR[1], pR[2], pR[3]);
        printf("            (%.12lg %.12lg %.12lg %.12lg)\n", 
                                pR[4], pR[5], pR[6], pR[7]);
        printf("e = (%.3lf %.3lf %.3lf %.3lf)\n", e[0],e[1], e[2], e[3]);
        printf("    (%.3lf %.3lf %.3lf %.3lf)\n", e[4],e[5], e[6], e[7]);
        printf("    (%.3lf %.3lf %.3lf %.3lf)\n", e[8],e[9], e[10], e[11]);
        printf("    (%.3lf %.3lf %.3lf %.3lf)\n", e[12],e[13], e[14], e[15]);
        printf("ie = (%.3lf %.3lf %.3lf %.3lf)\n", ie[0],ie[1], ie[2], ie[3]);
        printf("    (%.3lf %.3lf %.3lf %.3lf)\n", ie[4],ie[5], ie[6], ie[7]);
        printf("    (%.3lf %.3lf %.3lf %.3lf)\n", ie[8],ie[9], ie[10], ie[11]);
        printf("    (%.3lf %.3lf %.3lf %.3lf)\n", ie[12],ie[13], ie[14], ie[15]);
        printf("Uhll:  %.10lg %.10lg %.10lg %.10lg\n",
                Uhll[0], Uhll[1], Uhll[2], Uhll[3]);
        printf("Uhlld: %.10lg %.10lg %.10lg %.10lg\n",
                U[0], U[1], U[2], U[3]);
        printf("Uhll:  %.10lg %.10lg %.10lg %.10lg\n",
                Uhll[4], Uhll[5], Uhll[6], Uhll[7]);
        printf("Uhlld: %.10lg %.10lg %.10lg %.10lg\n",
                U[4], U[5], U[6], U[7]);
        printf("Fhll:  %.10lg %.10lg %.10lg %.10lg\n",
                Fhll[0], Fhll[1], Fhll[2], Fhll[3]);
        printf("Fhlld: %.10lg %.10lg %.10lg %.10lg\n",
                F[0], F[1], F[2], F[3]);
        printf("Fhll:  %.10lg %.10lg %.10lg %.10lg\n",
                Fhll[4], Fhll[5], Fhll[6], Fhll[7]);
        printf("Fhlld: %.10lg %.10lg %.10lg %.10lg\n",
                F[4], F[5], F[6], F[7]);

    }
    if(DEBUG2)
        printf("Done get_Ustar_HLLD.\n");
}

void calc_tetrads(double al, double *be, double *gam, double *igam, 
                    double sqrtgam, int xx, double *e, double *ie)
{
    if(DEBUG2)
        printf("  In calc_tetrads\n");
    //Here sqrtgam is the ACTUAL sqrtgam (not the sqrtgam/r used elsewhere).
    int yy = (xx+1)%3;
    int zz = (xx+2)%3;

    double ial = 1.0/al;
    double xfac = 1.0/sqrt(igam[3*xx+xx]); // 1.0 / sqrt(gam^{11})

    // e[4*i+mu] = e_(\hat{i})^\mu
    // ie[4*i+mu] = e^(\hat{i})_\mu

    //e_(t)
    e[4*0+0] = ial;
    e[4*0+xx+1] = -be[xx]*ial;
    e[4*0+yy+1] = -be[yy]*ial;
    e[4*0+zz+1] = -be[zz]*ial;
    //e_(x)
    e[4*1+0] = 0.0;
    e[4*1+xx+1] = igam[3*xx+xx] * xfac;
    e[4*1+yy+1] = igam[3*xx+yy] * xfac;
    e[4*1+zz+1] = igam[3*xx+zz] * xfac;
    //e_(y)
    e[4*2+0] = 0.0;
    e[4*2+xx+1] = 0.0;
    e[4*2+yy+1] = sqrt(gam[3*zz+zz]/(igam[3*xx+xx]*sqrtgam*sqrtgam));
    e[4*2+zz+1] = -gam[3*yy+zz] / (sqrtgam*sqrt(igam[3*xx+xx]*gam[3*zz+zz]));
    //e_(z)
    e[4*3+0] = 0.0;
    e[4*3+xx+1] = 0.0;
    e[4*3+yy+1] = 0.0;
    e[4*3+zz+1] = 1.0/sqrt(gam[3*zz+zz]);

    //e^(t)
    ie[4*0+0] = al;
    ie[4*0+xx+1] = 0.0;
    ie[4*0+yy+1] = 0.0;
    ie[4*0+zz+1] = 0.0;
    //e^(x)
    ie[4*1+0] = 0.0;
    ie[4*1+xx+1] = gam[3*xx+xx]*e[4*1+xx+1] + gam[3*xx+yy]*e[4*1+yy+1]
                    + gam[3*xx+zz]*e[4*1+zz+1];
    ie[4*1+yy+1] = gam[3*yy+xx]*e[4*1+xx+1] + gam[3*yy+yy]*e[4*1+yy+1]
                    + gam[3*yy+zz]*e[4*1+zz+1];
    ie[4*1+zz+1] = gam[3*zz+xx]*e[4*1+xx] + gam[3*zz+yy]*e[4*1+yy]
                    + gam[3*zz+zz]*e[4*1+zz];
    //e^(y)
    ie[4*2+0] = 0.0;
    ie[4*2+xx+1] = 0.0;
    ie[4*2+yy+1] = gam[3*yy+yy]*e[4*2+yy+1] + gam[3*yy+zz]*e[4*2+zz+1];
    ie[4*2+zz+1] = gam[3*zz+yy]*e[4*2+yy+1] + gam[3*zz+zz]*e[4*2+zz+1];
    //e^(z)
    ie[4*3+0] = 0.0;
    ie[4*3+xx+1] = 0.0;
    ie[4*3+yy+1] = 0.0;
    ie[4*3+zz+1] = gam[3*zz+zz]*e[4*3+zz+1];
    
    if(DEBUG2)
        printf("  Done calc_tetrads\n");
}

void calc_UF(double rho, double P, double *u, double *B, double *U, double *F)
{
    //Calculates U&F suitable for relativistic HLLD calculation,
    //Assumes u and B are in orthonormal frame.

    if(DEBUG2)
        printf("  In calc_UF\n");

    double rhoh = rho + gamma_law/(gamma_law-1)*P;
    double b[4], b2;
    b[0] = u[1]*B[0]+u[2]*B[1]+u[3]*B[2];
    b[1] = (B[0] + b[0]*u[1]) / u[0];
    b[2] = (B[1] + b[0]*u[2]) / u[0];
    b[3] = (B[2] + b[0]*u[3]) / u[0];
    b2 = -b[0]*b[0] + b[1]*b[1] + b[2]*b[2] + b[3]*b[3];
    U[DD] = rho*u[0];
    U[EN] = (rhoh+b2)*u[0]*u[0] - (P+0.5*b2) - b[0]*b[0];
    U[MX] = (rhoh+b2)*u[0]*u[1] - b[0]*b[1];
    U[MY] = (rhoh+b2)*u[0]*u[2] - b[0]*b[2];
    U[MZ] = (rhoh+b2)*u[0]*u[3] - b[0]*b[3];
    U[BX] = B[0];
    U[BY] = B[1];
    U[BZ] = B[2];

    F[DD] = rho*u[1];
    F[EN] = U[MX];
    F[MX] = (rhoh+b2)*u[1]*u[1] + (P+0.5*b2) - b[1]*b[1];
    F[MY] = (rhoh+b2)*u[2]*u[1] - b[2]*b[1];
    F[MZ] = (rhoh+b2)*u[3]*u[1] - b[3]*b[1];
    F[BX] = 0.0;
    F[BY] = (B[1]*u[1] - B[0]*u[2]) / u[0];
    F[BZ] = (B[2]*u[1] - B[0]*u[3]) / u[0];

    if(DEBUG2)
        printf("  Done calc_UF\n");
}

int solve_HLLD_SR(double sL, double sR, double Bx, double *UL, double *FL,
                    double *UR, double *FR, double w, double *pL, double *pR, 
                    double *U, double *F)
{
    if(DEBUG2)
        printf("  In solve_HLLD_SR\n");

    double Ps;

    double RL[NUM_C], RR[NUM_C], Uhll[NUM_C], Fhll[NUM_C];
    int q;
    for(q = 0; q<NUM_C; q++)
    {
        RL[q] = sL*UL[q] - FL[q];
        RR[q] = sR*UR[q] - FR[q];
        Uhll[q] = (sR*UR[q] - sL*UL[q] + FL[q]-FR[q]) / (sR - sL);
        Fhll[q] = (sR*FL[q] - sL*FR[q] + sL*sR*(UR[q]-UL[q])) / (sR - sL);
    }

    double prim[5];
    for(q=0; q<5; q++)
        prim[q] = 0.5*(pL[q]+pR[q]);
    double Pshll = cons2prim_hlld(Uhll, prim);
    
    //Eq. (55) of MUB
    double b = Uhll[EN]-Fhll[MX];
    double c = Uhll[MX]*Fhll[EN] - Fhll[MX]*Uhll[EN];
    double Ps0 = 0.5*(-b + sqrt(b*b-4*c));
    if(Bx*Bx < 0.1*Pshll)
        Ps = Ps0;
    else
        Ps = Pshll;

    if(DEBUG2 || DEBUG4)
        printf("    Ps0 = %.12lg (PsHLL=%.16lg PsB0=%.16lg)\n", Ps,Pshll,Ps0); 

    int err, mag;
    mag = Bx*Bx > 1.0e-16 * Ps ? 1 : 0;
    
    double BaL[3], Bc[3], BaR[3];
    double vaL[3], vcL[3], vcR[3], vaR[3];
    double KL[3], KR[3];
    double enthL, enthR;
    double saL, sc, saR;

    //TODO: THIS IS SUCH A HACK. calc_Ps should be able to handle B=0.
    int success;
    if(mag)
    {
        err = calc_Ps(RL, RR, Bx, sL, sR, &Ps);
        calc_state(Ps, RL, RR, Bx, sL, sR, vaL, vaR, BaL, BaR, &enthL, &enthR, 
                KL, KR, Bc, vcL, vcR);
        success = checkSolution(Ps, sL, sR, vaL, vaR, BaL, BaR, enthL, enthR, 
                                KL, KR, Bc, vcL, vcR);
        if(!success)
        {
            if(DEBUG4)
                printf("Second Attempt.\n");
            if(Bx*Bx < 0.1*Pshll)
                Ps = Pshll;
            else
                Ps = Ps0;
            err = calc_Ps(RL, RR, Bx, sL, sR, &Ps);
            calc_state(Ps, RL, RR, Bx, sL, sR, vaL, vaR, BaL, BaR, 
                        &enthL, &enthR, KL, KR, Bc, vcL, vcR);
            success = checkSolution(Ps, sL, sR, vaL, vaR, BaL, BaR, enthL, 
                                    enthR, KL, KR, Bc, vcL, vcR);
        }
        if(!success)
        {
            if(DEBUG4)
                printf("Third Attempt.\n");

            double Ap = RL[MX]-sL*RL[EN];
            double C = RL[MY]*RL[BY] + RL[MZ]*RL[BZ];
            double G = RL[BY]*RL[BY] + RL[BZ]*RL[BZ];
            double a = -sL*(1.0-sL*sL);
            double b = (sL*Bx*Bx - RL[EN])*(1.0-sL*sL) - sL*(Ap+G);
            double c = Bx*(sL*Bx*Ap + C) - (Ap+G)*RL[EN]; 

            double desc = b*b-4*a*c;
            if(desc < 0)
                success = 0;
            else
            {
                double Psing1 = (-b - sqrt(b*b-4*a*c)) / (2*a);
                double Psing2 = (-b + sqrt(b*b-4*a*c)) / (2*a);
                double del = 1.0e-1;
                if(Psing1 > 0.0)
                {
                    if(DEBUG4)
                        printf("PL1-.\n");
                    Ps = Psing1*(1.0+del);
                    err = calc_Ps(RL, RR, Bx, sL, sR, &Ps);
                    calc_state(Ps, RL, RR, Bx, sL, sR, vaL, vaR, BaL, BaR, 
                                &enthL, &enthR, KL, KR, Bc, vcL, vcR);
                    success = checkSolution(Ps, sL, sR, vaL, vaR, BaL, BaR, 
                                            enthL, enthR, KL, KR, Bc, 
                                            vcL, vcR);
                    if(!success)
                    {
                        if(DEBUG4)
                            printf("PL1+.\n");
                        Ps = Psing1*(1.0-del);
                        err = calc_Ps(RL, RR, Bx, sL, sR, &Ps);
                        calc_state(Ps, RL, RR, Bx, sL, sR, vaL, vaR, BaL, BaR, 
                                    &enthL, &enthR, KL, KR, Bc, vcL, vcR);
                        success = checkSolution(Ps, sL, sR, vaL, vaR, BaL, 
                                                BaR, 
                                                enthL, enthR, KL, KR, Bc, 
                                                vcL, vcR);
                    }
                }
                if(!success && Psing2 > 0.0)
                {
                    if(DEBUG4)
                        printf("PL2-.\n");
                    Ps = Psing2*(1.0+del);
                    err = calc_Ps(RL, RR, Bx, sL, sR, &Ps);
                    calc_state(Ps, RL, RR, Bx, sL, sR, vaL, vaR, BaL, BaR, 
                                &enthL, &enthR, KL, KR, Bc, vcL, vcR);
                    success = checkSolution(Ps, sL, sR, vaL, vaR, BaL, BaR, 
                                            enthL, enthR, KL, KR, Bc, 
                                            vcL, vcR);
                    if(!success)
                    {
                        if(DEBUG4)
                            printf("PL2+.\n");
                        Ps = Psing2*(1.0-del);
                        err = calc_Ps(RL, RR, Bx, sL, sR, &Ps);
                        calc_state(Ps, RL, RR, Bx, sL, sR, vaL, vaR, BaL, BaR, 
                                    &enthL, &enthR, KL, KR, Bc, vcL, vcR);
                        success = checkSolution(Ps, sL, sR, vaL, vaR, BaL, 
                                                BaR, 
                                                enthL, enthR, KL, KR, Bc, 
                                                vcL, vcR);
                    }
                }
            }
        }
        if(!success)
        {
            if(DEBUG4)
                printf("Fourth Attempt.\n");

            double Ap = RR[MX]-sR*RR[EN];
            double C = RR[MY]*RR[BY] + RR[MZ]*RR[BZ];
            double G = RR[BY]*RR[BY] + RR[BZ]*RR[BZ];
            double a = -sR*(1.0-sR*sR);
            double b = (sR*Bx*Bx - RR[EN])*(1.0-sR*sR) - sR*(Ap+G);
            double c = Bx*(sR*Bx*Ap + C) - (Ap+G)*RR[EN]; 

            double desc = b*b-4*a*c;
            if(desc < 0)
                success = 0;
            else
            {
                double Psing1 = (-b + sqrt(b*b-4*a*c)) / (2*a);
                double Psing2 = (-b - sqrt(b*b-4*a*c)) / (2*a);
                double del = 1.0e-1;
                if(Psing1 > 0.0)
                {
                    if(DEBUG4)
                        printf("PR1-.\n");
                    Ps = Psing1*(1.0+del);
                    err = calc_Ps(RL, RR, Bx, sL, sR, &Ps);
                    calc_state(Ps, RL, RR, Bx, sL, sR, vaL, vaR, BaL, BaR, 
                                &enthL, &enthR, KL, KR, Bc, vcL, vcR);
                    success = checkSolution(Ps, sL, sR, vaL, vaR, BaL, BaR, 
                                            enthL, enthR, KL, KR, Bc, 
                                            vcL, vcR);
                    if(!success)
                    {
                        if(DEBUG4)
                            printf("PR1+.\n");
                        Ps = Psing1*(1.0-del);
                        err = calc_Ps(RL, RR, Bx, sL, sR, &Ps);
                        calc_state(Ps, RL, RR, Bx, sL, sR, vaL, vaR, BaL, BaR, 
                                    &enthL, &enthR, KL, KR, Bc, vcL, vcR);
                        success = checkSolution(Ps, sL, sR, vaL, vaR, BaL, 
                                                BaR, 
                                                enthL, enthR, KL, KR, Bc, 
                                                vcL, vcR);
                    }
                }
                if(!success && Psing2 > 0.0)
                {
                    if(DEBUG4)
                        printf("PR2-.\n");
                    Ps = Psing2*(1.0+del);
                    err = calc_Ps(RL, RR, Bx, sL, sR, &Ps);
                    calc_state(Ps, RL, RR, Bx, sL, sR, vaL, vaR, BaL, BaR, 
                                &enthL, &enthR, KL, KR, Bc, vcL, vcR);
                    success = checkSolution(Ps, sL, sR, vaL, vaR, BaL, BaR, 
                                            enthL, enthR, KL, KR, Bc, 
                                            vcL, vcR);
                    if(!success)
                    {
                        if(DEBUG4)
                            printf("PR2+.\n");
                        Ps = Psing2*(1.0-del);
                        err = calc_Ps(RL, RR, Bx, sL, sR, &Ps);
                        calc_state(Ps, RL, RR, Bx, sL, sR, vaL, vaR, BaL, BaR, 
                                    &enthL, &enthR, KL, KR, Bc, vcL, vcR);
                        success = checkSolution(Ps, sL, sR, vaL, vaR, BaL, 
                                                BaR, 
                                                enthL, enthR, KL, KR, Bc, 
                                                vcL, vcR);
                    }
                }
            }
        }
    }
    else
    {
        calc_state(Ps, RL, RR, Bx, sL, sR, vaL, vaR, BaL, BaR, 
                    &enthL, &enthR, KL, KR, Bc, vcL, vcR);
        success = checkSolution(Ps, sL, sR, vaL, vaR, BaL, BaR, enthL, 
                                enthR, KL, KR, Bc, vcL, vcR);
    }
    if(DEBUG2 || DEBUG4)
        printf("    Ps = %.12lg\n", Ps); 



    saL = KL[0];
    saR = KR[0];
    sc = 0.5*(vcL[0] + vcR[0]);

    if(DEBUG4)
    {
        printf("UL: %.16lg %.16lg %.16lg %.16lg\n", 
                UL[0], UL[1], UL[2], UL[3]);
        printf("    %.16lg %.16lg %.16lg %.16lg\n", 
                UL[4], UL[5], UL[6], UL[7]);
        printf("UR: %.16lg %.16lg %.16lg %.16lg\n", 
                UR[0], UR[1], UR[2], UR[3]);
        printf("    %.16lg %.16lg %.16lg %.16lg\n", 
                UR[4], UR[5], UR[6], UR[7]);
        printf("FL: %.16lg %.16lg %.16lg %.16lg\n", 
                FL[0], FL[1], FL[2], FL[3]);
        printf("    %.16lg %.16lg %.16lg %.16lg\n", 
                FL[4], FL[5], FL[6], FL[7]);
        printf("FR: %.16lg %.16lg %.16lg %.16lg\n", 
                FR[0], FR[1], FR[2], FR[3]);
        printf("    %.16lg %.16lg %.16lg %.16lg\n", 
                FR[4], FR[5], FR[6], FR[7]);
        printf("Bx = %.16lg, sL = %.16lg, sR = %.16lg\n", Bx, sL, sR);
        printf("enth: %.12lg %.12lg\n", enthL, enthR);

    printf("B: (%.6lg %.6lg %.6lg) (%.6lg %.6lg %.6lg) (%.6lg %.6lg %.6lg)\n",
            BaL[0], BaL[1], BaL[2], Bc[0], Bc[1], Bc[2], 
            BaR[0], BaR[1], BaR[2]);
    printf("K: (%.12lg %.6lg %.6lg) (%.12lg %.6lg %.6lg)\n",
            KL[0], KL[1], KL[2], KR[0], KR[1], KR[2]);
    printf("va: (%.6lg %.6lg %.6lg) (%.6lg %.6lg %.6lg)\n",
            vaL[0], vaL[1], vaL[2], vaR[0], vaR[1], vaR[2]);
    printf("vc: (%.6lg %.6lg %.6lg) (%.6lg %.6lg %.6lg)\n",
            vcL[0], vcL[1], vcL[2], vcR[0], vcR[1], vcR[2]);
    }

    int zone = -1;

    if(success)
    {
        if(w < sL)
        {
            if(DEBUG3 || DEBUG4)
                printf("    State  L!! (w=%.8lg)\n", w);
            for(q=0; q<NUM_C; q++)
            {
                U[q] = UL[q];
                F[q] = FL[q];
            }
            zone = 0;
        }
        else if (w < saL)
        {
            if(DEBUG3 || DEBUG4)
                printf("    State aL!! (w=%.8lg)\n", w);
            calc_Ua(Ps, sL, Bx, RL, vaL, BaL, U);
            for(q=0; q<NUM_C; q++)
                F[q] = FL[q] + sL*(U[q]-UL[q]);
            zone = 1;
        }
        else if(w < sc && mag)
        {
            if(DEBUG3 || DEBUG4)
                printf("    State cL!! (w=%.8lg)\n", w);
            double UaL[NUM_C];
            calc_Ua(Ps, sL, Bx, RL, vaL, BaL, UaL);
            calc_Uc(Ps, saL, Bx, vaL, vcL, Bc, UaL, U);
            for(q=0; q<NUM_C; q++)
                F[q] = FL[q] + sL*(UaL[q]-UL[q]) + saL*(U[q]-UaL[q]);
            zone = 2;
        }
        else if(w < saR && mag)
        {
            if(DEBUG3 || DEBUG4)
                printf("    State cR!! (w=%.8lg)\n", w);
            double UaR[NUM_C];
            calc_Ua(Ps, sR, Bx, RR, vaR, BaR, UaR);
            calc_Uc(Ps, saR, Bx, vaR, vcR, Bc, UaR, U);
            for(q=0; q<NUM_C; q++)
                F[q] = FR[q] + sR*(UaR[q]-UR[q]) + saR*(U[q]-UaR[q]);
            zone = 3;
        }
        else if(w < sR)
        {
            if(DEBUG3 || DEBUG4)
                printf("    State aR!! (w=%.8lg)\n", w);
            calc_Ua(Ps, sR, Bx, RR, vaR, BaR, U);
            for(q=0; q<NUM_C; q++)
                F[q] = FR[q] + sR*(U[q]-UR[q]);
            zone = 4;
        }
        else
        {
            if(DEBUG3 || DEBUG4)
                printf("    State  R!! (w=%.8lg)\n", w);
            for(q=0; q<NUM_C; q++)
            {
                U[q] = UR[q];
                F[q] = FR[q];
            }
            zone = 5;
        }
    }
    else
    {
        if(w < sL)
        {
            if(DEBUG3 || DEBUG4)
                printf("    State  L!! (w=%.8lg)\n", w);
            for(q=0; q<NUM_C; q++)
            {
                U[q] = UL[q];
                F[q] = FL[q];
            }
            zone = -3;
        }
        else if (w < sR)
        {
            if(DEBUG3 || DEBUG4)
                printf("    State HLL!! (w=%.8lg)\n", w);
            for(q=0; q<NUM_C; q++)
            {
                U[q] = Uhll[q];
                F[q] = Fhll[q];
            }
            zone = -2;
        }
        else
        {
            if(DEBUG3 || DEBUG4)
                printf("    State  R!! (w=%.8lg)\n", w);
            for(q=0; q<NUM_C; q++)
            {
                U[q] = UR[q];
                F[q] = FR[q];
            }
            zone = -1;
        }
    }

    if(DEBUG4)
        printf("Wavespeeds: %.8lg %.8lg %.8lg %.8lg %.8lg\n",
                            sL, saL, sc, saR, sR);

    if(DEBUG2)
        printf("  Done solve_HLLD_SR\n");

    return zone;
}

int checkSolution(double Ps, double sL, double sR, double *vaL, double *vaR, 
                    double *BaL, double *BaR, double enthL, double enthR, 
                    double *KaL, double *KaR, double *Bc, double *vcL, 
                    double *vcR)
{
    int success = 1;

    if(Ps != Ps)
        success = 0;
    if(enthL < Ps || enthR < Ps)
        success = 0;
    if(enthL < Ps || enthR < Ps)
        success = 0;
    if(fabs(vaL[0]) > 1.0 || fabs(vaL[1]) > 1.0 || fabs(vaL[2]) > 1.0)
        success = 0;
    if(fabs(vaR[0]) > 1.0 || fabs(vaR[1]) > 1.0 || fabs(vaR[2]) > 1.0)
        success = 0;
    if(vaL[0] < sL || vaR[0] > sR)
        success = 0;
    if(fabs(vcL[0]) > 1.0 || fabs(vcL[1]) > 1.0 || fabs(vcL[2]) > 1.0)
        success = 0;
    if(fabs(vcR[0]) > 1.0 || fabs(vcR[1]) > 1.0 || fabs(vcR[2]) > 1.0)
        success = 0;
    if(vcL[0] < sL || vcR[0] > sR)
        success = 0;
    if(KaL[0] < sL || KaR[0] > sR)
        success = 0;

    return success;
}

void calc_va(double p, double s, double Bx, double *R, double *va, double *dva)
{
    //Eq. (23-25) of MUB
    if(DEBUG3)
        printf("      In calc_va.\n");
    double A = R[MX] - s*R[EN] + p*(1.0-s*s);
    double G = R[BY]*R[BY] + R[BZ]*R[BZ];
    double C = R[MY]*R[BY] + R[MZ]*R[BZ];
    double Q = -A - G + Bx*Bx*(1.0-s*s);
    double iX = 1.0/(Bx * (s*A*Bx+C) - (A+G) * (s*p+R[EN]));

    va[0] = iX * (Bx * (A*Bx+s*C) - (A+G) * (p+R[MX]));
    va[1] = iX * (Q*R[MY] + R[BY] * (C+Bx*(s*R[MX]-R[EN])));
    va[2] = iX * (Q*R[MZ] + R[BZ] * (C+Bx*(s*R[MX]-R[EN])));

    if(dva != NULL)
    {
        double dA = 1.0-s*s;
        double dQ = -dA;
        double diX = -iX*iX*(s*Bx*Bx*dA - s*(A+G) - dA*(s*p+R[EN]));
        dva[0] = diX*(Bx * (A*Bx+s*C) - (A+G) * (p+R[MX]))
                    + iX*(Bx*Bx*dA - dA*(p+R[MX]) - (A+G));
        dva[1] = diX * (Q*R[MY] + R[BY] * (C+Bx*(s*R[MX]-R[EN])))
                    + iX*(dQ*R[MY]);
        dva[2] = diX * (Q*R[MZ] + R[BZ] * (C+Bx*(s*R[MX]-R[EN])))
                    + iX*(dQ*R[MZ]);
    }
    if(DEBUG3)
        printf("      Done calc_va.\n");
}

void calc_Ba(double s, double Bx, double *R, double *va, double *Ba, 
                double *dva, double *dBa)
{
    //Eq. (21) of MUB
    if(DEBUG3)
        printf("      In calc_Ba.\n");
    Ba[0] = Bx;
    Ba[1] = (R[BY] - Bx*va[1]) / (s - va[0]);
    Ba[2] = (R[BZ] - Bx*va[2]) / (s - va[0]);

    if(dBa != NULL)
    {
        double ids2 = 1.0/((s-va[0])*(s-va[0]));
        dBa[0] = 0.0;
        dBa[1] = (-Bx*dva[1]*(s-va[0]) + (R[BY]-Bx*va[1])*dva[0]) * ids2;
        dBa[2] = (-Bx*dva[2]*(s-va[0]) + (R[BZ]-Bx*va[2])*dva[0]) * ids2;
    }
    if(DEBUG3)
        printf("      Done calc_Ba.\n");
}

void calc_entha(double p, double s, double *R, double *va, double *enth,
                    double *dva, double *denth)
{
    //Eq. (31) of MUB
    if(DEBUG3)
        printf("      In calc_entha.\n");
    *enth =  p + (R[EN] - va[0]*R[MX]-va[1]*R[MY]-va[2]*R[MZ]) / (s - va[0]);

    if(denth != NULL)
        *denth = 1.0 + ((-dva[0]*R[MX]-dva[1]*R[MY]-dva[2]*R[MZ])*(s-va[0])
                        + (R[EN]-va[0]*R[MX]-va[1]*R[MY]-va[2]*R[MZ])*dva[0])
                        / ((s-va[0])*(s-va[0]));
    if(DEBUG3)
        printf("      Done calc_entha.\n");
}

void calc_Ka(double p, double s, double Bx, double enth, double *R, double *Ka,
                double denth, double *dKa, int LR)
{
    //Eq. (43) of MUB
    //LR == 1 for Right, LR == -1 for Left
    if(DEBUG3)
        printf("      In calc_Ka.\n");
    
    int sgn = Bx > 0 ? LR : -LR;
    double eta = sgn*sqrt(enth);
    double denom = 1.0 / (s*p + R[EN] + Bx*eta);

    Ka[0] = (R[MX] + p + R[BX]*eta) * denom;
    Ka[1] = (R[MY]     + R[BY]*eta) * denom;
    Ka[2] = (R[MZ]     + R[BZ]*eta) * denom;

    if(dKa != NULL)
    {
        double deta = 0.5/eta * denth;
        double ddenom = -denom*denom*(s + Bx*deta);

        dKa[0] = (1.0 + R[BX]*deta)*denom + (R[MX]+p+R[BX]*eta)*ddenom;
        dKa[1] = (R[BY]*deta)*denom + (R[MY]+R[BY]*eta)*ddenom;
        dKa[2] = (R[BZ]*deta)*denom + (R[MZ]+R[BZ]*eta)*ddenom;
    } 
    if(DEBUG3)
        printf("      Done calc_Ka.\n");
}

void cont_func(double Bx, double *BaL, double *BaR,
                double *vaL, double *vaR, double *KaL, double *KaR,
                double enthL, double enthR, double *f,
                double *dBaL, double *dBaR, double *dvaL, double *dvaR,
                double *dKaL, double *dKaR, double denthL, double denthR,
                double *df)
{
    //Eq. (48) of MUB
    if(DEBUG3)
        printf("      In cont_func.\n");

    double saL = KaL[0];
    double saR = KaR[0];
    double dsaL = dKaL[0];
    double dsaR = dKaR[0];
    
    // Bch = \hat{B_c}, numerator of Eq. 45 for Bc. Recalculating here is more
    //  stable if B~0.
    double Bch[3];
    Bch[0] = (saR-saL)*Bx;
    Bch[1] =  (BaR[1]*(saR-vaR[0]) + Bx*vaR[1]) 
            - (BaL[1]*(saL-vaL[0]) + Bx*vaL[1]);
    Bch[2] =  (BaR[2]*(saR-vaR[0]) + Bx*vaR[2]) 
            - (BaL[2]*(saL-vaL[0]) + Bx*vaL[2]);

    double KaL2 = KaL[0]*KaL[0] + KaL[1]*KaL[1] + KaL[2]*KaL[2];
    double KaR2 = KaR[0]*KaR[0] + KaR[1]*KaR[1] + KaR[2]*KaR[2];
    double KBL = KaL[0]*Bch[0] + KaL[1]*Bch[1] + KaL[2]*Bch[2];
    double KBR = KaR[0]*Bch[0] + KaR[1]*Bch[1] + KaR[2]*Bch[2];

    int sgn = Bx>0 ? 1 : -1;
    double etaL = -sgn*sqrt(enthL);
    double etaR =  sgn*sqrt(enthR);

    double delKx = (KaR[0]-KaL[0]);

    // Eq. 49
    double YL = (1.0 - KaL2) / (etaL*delKx - KBL);
    double YR = (1.0 - KaR2) / (etaR*delKx - KBR);

    //Eq. 48
    *f = (KaR[0]-KaL[0]) * (1.0 - Bx*(YR-YL));

    if(df != NULL)
    {
        double dBch[3];
        dBch[0] = (dsaR - dsaL) * Bx;
        dBch[1] = (dBaR[1]*(saR-vaR[0]) + BaR[1]*(dsaR-dvaR[0]) + Bx*dvaR[1])
                - (dBaL[1]*(saL-vaL[0]) + BaL[1]*(dsaL-dvaL[0]) + Bx*dvaL[1]);
        dBch[2] = (dBaR[2]*(saR-vaR[0]) + BaR[2]*(dsaR-dvaR[0]) + Bx*dvaR[2])
                - (dBaL[2]*(saL-vaL[0]) + BaL[2]*(dsaL-dvaL[0]) + Bx*dvaL[2]);
        
        double dKaL2 = 2*(KaL[0]*dKaL[0] + KaL[1]*dKaL[1] + KaL[2]*dKaL[2]);
        double dKaR2 = 2*(KaR[0]*dKaR[0] + KaR[1]*dKaR[1] + KaR[2]*dKaR[2]);
        double dKBL = dKaL[0]*Bch[0] + KaL[0]*dBch[0]
                    + dKaL[1]*Bch[1] + KaL[1]*dBch[1]
                    + dKaL[2]*Bch[2] + KaL[2]*dBch[2];
        double dKBR = dKaR[0]*Bch[0] + KaR[0]*dBch[0]
                    + dKaR[1]*Bch[1] + KaR[1]*dBch[1]
                    + dKaR[2]*Bch[2] + KaR[2]*dBch[2];

        double detaL = 0.5/etaL * denthL;
        double detaR = 0.5/etaR * denthR;

        double ddelKx = dKaR[0] - dKaL[0];
        
        double dYL = (-dKaL2*(etaL*delKx-KBL)
                        - (1-KaL2)*(detaL*delKx + etaL*ddelKx - dKBL))
                    / ((etaL*delKx-KBL)*(etaL*delKx-KBL));
        double dYR = (-dKaR2*(etaR*delKx-KBR)
                        - (1-KaR2)*(detaR*delKx + etaR*ddelKx - dKBR))
                    / ((etaR*delKx-KBR)*(etaR*delKx-KBR));

        *df = ddelKx * (1.0-Bx*(YR-YL)) + delKx * (-Bx*(dYR-dYL));
    }
    if(DEBUG3)
        printf("      Done cont_func.\n");
}

void calc_Bc(double Bx, double *BaL, double *BaR, double *vaL, double *vaR, 
                double saL, double saR, double *Bc)
{
    //Eq. (45) of MUB
    
    if(DEBUG2)
        printf("    In calc_Bc.\n");

    Bc[0] = Bx;
    Bc[1] = ((BaR[1]*(saR-vaR[0])+Bx*vaR[1])
                - (BaL[1]*(saL-vaL[0]) + Bx*vaL[1]))
            / (saR - saL);
    Bc[2] = ((BaR[2]*(saR-vaR[0])+Bx*vaR[2])
                - (BaL[2]*(saL-vaL[0]) + Bx*vaL[2]))
            / (saR - saL);

    if(Bc[1] != Bc[1])
        Bc[1] = 0.0;
    if(Bc[2] != Bc[2])
        Bc[2] = 0.0;
    
    if(DEBUG2)
        printf("    Done calc_Bc.\n");
}

void calc_vc(double Bx, double *Ka, double *Bc, double enth, double *vc,
                int LR)
{
    int sgn = Bx > 0 ? LR : -LR;
    double eta = sgn*sqrt(enth);

    double K2 = Ka[0]*Ka[0] + Ka[1]*Ka[1] + Ka[2]*Ka[2];
    double KB = Ka[0]*Bx + Ka[1]*Bc[1] + Ka[2]*Bc[2];

    vc[0] = Ka[0] - Bx * (1.0-K2) / (eta - KB);
    vc[1] = Ka[1] - Bc[1] * (1.0-K2) / (eta - KB);
    vc[2] = Ka[2] - Bc[2] * (1.0-K2) / (eta - KB);
}

int calc_Ps(double *RL, double *RR, double Bx, double sL, double sR, 
                double *Ps)
{ 
    if(DEBUG2)
        printf("    In calc_Ps\n");

    int i;
    double p0, p, p1, dp;
    double f, df;

    p0 = *Ps;

    p1 = p0;

    i = 0;
    do{
        p = p1;
        double vaL[3], vaR[3], dvaL[3], dvaR[3];
        double BaL[3], BaR[3], dBaL[3], dBaR[3];
        double enthL, enthR, denthL, denthR;
        double KaL[3], KaR[3], dKaL[3], dKaR[3];

        calc_va(p, sL, Bx, RL, vaL, dvaL);
        calc_va(p, sR, Bx, RR, vaR, dvaR);
        calc_Ba(sL, Bx, RL, vaL, BaL, dvaL, dBaL);
        calc_Ba(sR, Bx, RR, vaR, BaR, dvaR, dBaR);
        calc_entha(p, sL, RL, vaL, &enthL, dvaL, &denthL);
        calc_entha(p, sR, RR, vaR, &enthR, dvaR, &denthR);
        calc_Ka(p, sL, Bx, enthL, RL, KaL, denthL, dKaL, -1);
        calc_Ka(p, sR, Bx, enthR, RR, KaR, denthR, dKaR, +1);

        cont_func(Bx, BaL, BaR, vaL, vaR, KaL, KaR, enthL, enthR, &f,
                    dBaL, dBaR, dvaL, dvaR, dKaL, dKaR, denthL, denthR, &df);

        dp = -f/df;

        p1 = p + dp;

        if(p1 < 0.0)
            p1 = 0.5*p;

        if(enthL < 0.0 || enthR < 0.0)
            p1 = 0.8*p;

        if(DEBUG4)
        {
            printf("      %d: p=%.12lg dp=%.12lg f=%.12lg df=%.12lg\n",
                    i, p, dp, f, df);
            //printf("          dBaL=(%.6lg %.6lg %.6lg) dBaR=(%.6lg %.6lg %.6lg)\n", dBaL[0], dBaL[1], dBaL[2], dBaR[0], dBaR[1], dBaR[2]);
            //printf("          dvaL=(%.6lg %.6lg %.6lg) dvaR=(%.6lg %.6lg %.6lg)\n", dvaL[0], dvaL[1], dvaL[2], dvaR[0], dvaR[1], dvaR[2]);
            //printf("          dKaL=(%.6lg %.6lg %.6lg) dKaR=(%.6lg %.6lg %.6lg)\n", dKaL[0], dKaL[1], dKaL[2], dKaR[0], dKaR[1], dKaR[2]);
            //printf("          denthL=%.6lg denthR=%.6lg\n", denthL, denthR);
        }
        if(dp != dp && DEBUG)
        {
            if(!DEBUG4)
            {
                printf("      RL: %.16lg %.16lg %.16lg %.16lg\n", 
                        RL[0], RL[1], RL[2], RL[3]);
                printf("          %.16lg %.16lg %.16lg %.16lg\n", 
                        RL[4], RL[5], RL[6], RL[7]);
                printf("      RR: %.16lg %.16lg %.16lg %.16lg\n", 
                        RR[0], RR[1], RR[2], RR[3]);
                printf("          %.16lg %.16lg %.16lg %.16lg\n", 
                        RR[4], RR[5], RR[6], RR[7]);
                printf("      Bx = %.16lg, sL = %.16lg, sR = %.16lg\n",
                        Bx, sL, sR);
            }
            printf("        va: (%.6lg %.6lg %.6lg) (%.6lg %.6lg %.6lg)\n",
                    vaL[0], vaL[1], vaL[2], vaR[0], vaR[1], vaR[2]);
            printf("        Ba: (%.6lg %.6lg %.6lg) (%.6lg %.6lg %.6lg)\n",
                    BaL[0], BaL[1], BaL[2], BaR[0], BaR[1], BaR[2]);
            printf("        entha: %.12lg %.12lg\n", enthL, enthR);
            printf("        Ka: (%.6lg %.6lg %.6lg) (%.6lg %.6lg %.6lg)\n",
                    KaL[0], KaL[1], KaL[2], KaR[0], KaR[1], KaR[2]);
        }

        i++;
        
    }while(fabs(dp)/p > ACC && i < MAX_ITER);

    *Ps = p1;
    
    if(DEBUG2)
        printf("    Done calc_Ps.\n");

    return 0;
}

void calc_state(double p, double *RL, double *RR, double Bx, double sL, 
                    double sR, double *vaL, double *vaR, double *BaL, 
                    double *BaR, double *enthL, double *enthR, double *KaL,
                    double *KaR, double *Bc, double *vcL, double *vcR)
{
    if(DEBUG2)
        printf("    In calc_state.\n");
    calc_va(p, sL, Bx, RL, vaL, NULL);
    calc_va(p, sR, Bx, RR, vaR, NULL);
    calc_Ba(sL, Bx, RL, vaL, BaL, NULL, NULL);
    calc_Ba(sR, Bx, RR, vaR, BaR, NULL, NULL);
    calc_entha(p, sL, RL, vaL, enthL, NULL, NULL);
    calc_entha(p, sR, RR, vaR, enthR, NULL, NULL);
    calc_Ka(p, sL, Bx, *enthL, RL, KaL, 0.0, NULL, -1);
    calc_Ka(p, sR, Bx, *enthR, RR, KaR, 0.0, NULL, +1);
    calc_Bc(Bx, BaL, BaR, vaL, vaR, KaL[0], KaR[0], Bc);
    calc_vc(Bx, KaL, Bc, *enthL, vcL, -1);
    calc_vc(Bx, KaR, Bc, *enthR, vcR, +1);
    if(DEBUG2)
        printf("    Done calc_state.\n");
}

void calc_Ua(double p, double s, double Bx, double *R, double *va, double *Ba,
                double *Ua)
{
    //Eq. (32-34) of MUB
    double vB = va[0]*Bx+va[1]*Ba[1]+va[2]*Ba[2];

    Ua[DD] = R[DD] / (s - va[0]);
    Ua[EN] = (R[EN] + p*va[0] - vB*Bx) / (s-va[0]);
    Ua[MX] = (Ua[EN] + p) * va[0] - vB*Bx;
    Ua[MY] = (Ua[EN] + p) * va[1] - vB*Ba[1];
    Ua[MZ] = (Ua[EN] + p) * va[2] - vB*Ba[2];
    Ua[BX] = Bx;
    Ua[BY] = Ba[1];
    Ua[BZ] = Ba[2];
}

void calc_Uc(double p, double sa, double Bx, double *va, double *vc, 
                double *Bc, double *Ua, double *Uc)
{
    //Eq. (50-52) of MUB

    double vB = vc[0]*Bx + vc[1]*Bc[1] + vc[2]*Bc[2];

    Uc[DD] = Ua[DD] * (sa - va[0]) / (sa - vc[0]);
    Uc[EN] = (sa*Ua[EN] - Ua[MX] + p*vc[0] - vB*Bx) / (sa - vc[0]);
    Uc[MX] = (Uc[EN] + p)*vc[0] - vB*Bx;
    Uc[MY] = (Uc[EN] + p)*vc[1] - vB*Bc[1];
    Uc[MZ] = (Uc[EN] + p)*vc[2] - vB*Bc[2];
    Uc[BX] = Bx;
    Uc[BY] = Bc[1];
    Uc[BZ] = Bc[2];
}

double cons2prim_hlld(double *U, double *p)
{
    //if(isothermal)
    //    return cons2prim_isothermal_hlld(U, p);
        
    return cons2prim_noble2d_hlld(U, p);
}

double cons2prim_isothermal_hlld(double *U, double *p)
{
    double prec = 1.0e-10;
    double max_iter = 30;
    double Nextra = 2;

    double D = U[DD];
    double S[3] = {U[MX], U[MY], U[MZ]};
    double tau = U[EN];
    double B[3] = {U[BX], U[BY], U[BZ]};

    double cs2N = get_cs2(1.0);
    double P_o_rhoh = cs2N / gamma_law;
    double h = 1.0 / (1.0 - gamma_law * P_o_rhoh / (gamma_law-1.0));

    double s2 = (S[0]*S[0] + S[1]*S[1] + S[2]*S[2]) / (D*D*h*h);
    double B2 = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];
    double BS = B[0]*S[0] + B[1]*S[1] + B[2]*S[2];

    double Q = B2 / (D*h);
    double psi = BS*BS / (D*D*D*h*h*h);

    //Initial guess: previous wmo
    double u[3] = {p[URR], p[UPP], p[UZZ]};
    double u2 = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    double w = sqrt(1.0 + u2);
    double wmo0 = u2 / (w+1);

    //Run Newton-Raphson
    double wmo, wmo1;
    wmo1 = wmo0;
    
    int i = 0;
    int clean = -1;
    
    if(DEBUG4)
    {
        printf("H = %.12lg, s2 = %.12lg, Q = %.12lg, psi = %.12lg\n",
                    D*h, s2, Q, psi);
        printf("0: (%.12lg)\n", wmo1);
    }

    double c4 = 1.0;
    double c3 = 4 + 2*Q;
    double c2 = 5 - s2 + 6*Q + Q*Q;
    double c1 = 2 - 2*s2 + 4*Q + 2*Q*Q - 2*psi;
    double c0 = -s2 - 2*psi - Q*psi;

    double wmoMIN = 0.0;
    double wmoMAX = 1000.0;

    while(1)
    {
        wmo = wmo1;
        
        double f =  (((c4*wmo + c3)*wmo + c2)*wmo + c1)*wmo + c0;
        double df = ((4*c4*wmo + 3*c3)*wmo + 2*c2)*wmo + c1;

        wmo1  = wmo  - f/df;

        i++;

        if(wmo1 > wmoMAX || wmo1 < wmoMIN)
            wmo1 = 0.5*(wmoMAX+wmoMIN);
        else if (f > 0)
            wmoMAX = wmo;
        else if (f < 0)
            wmoMIN = wmo;

        double err = (wmo1-wmo) / (1.0+wmo);
        //if(err != err)
        //    printf("WHAT: v2=%.12lg, eta=%.12lg\n");

        if(DEBUG4)
        {
            printf("%d: (%.12lg) (%.12lg, %.12lg) %.12lg\n", 
                    i, wmo1, f, df, err);
        }

        if(fabs(err) < prec && clean < 0)
            clean = Nextra+1;
        if(clean >= 0)
            clean--;
        if(clean == 0 || i == max_iter)
            break;
    }

    if(i == max_iter && (DEBUG3 || DEBUG4) )
    {
        printf("ERROR: NR failed to converge in HLLD.  err = %.12lg\n", 
                fabs(wmo1-wmo)/(1+wmo));
        printf("    s2 = %.12lg, Q = %.12lg, psi = %.12lg\n", s2, Q, psi);
        printf("    wmo0 = %.12lg, wmo1 = %.12lg\n", wmo0, wmo1);
    }

    wmo = wmo1;

    //Prim recovery
    w = wmo+1.0;

    double rho = D / w;
    double Pp = P_o_rhoh * rho * h;
    
    u[0] = (S[0] + BS*B[0]/(D*h*w)) / (D*h + B2/w);
    u[1] = (S[1] + BS*B[1]/(D*h*w)) / (D*h + B2/w);
    u[2] = (S[2] + BS*B[2]/(D*h*w)) / (D*h + B2/w);

    double uB = u[0]*B[0] + u[1]*B[1] + u[2]*B[2];
    double b2 = (B2 + uB*uB) / (w*w);
    
    p[RHO] = rho;
    p[URR] = u[0];
    p[UPP] = u[1];
    p[UZZ] = u[2];
    p[PPP] = Pp;

    return Pp + 0.5*b2;
}

double cons2prim_noble2d_hlld(double *U, double *p)
{
    double prec = 1.0e-10;
    double max_iter = 30;
    double Nextra = 2;

    double D = U[DD];
    double S[3] = {U[MX], U[MY], U[MZ]};
    double tau = U[EN];
    double B[3] = {U[BX], U[BY], U[BZ]};

    double s2 = (S[0]*S[0] + S[1]*S[1] + S[2]*S[2]) / (D*D);
    double B2 = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];
    double BS = B[0]*S[0] + B[1]*S[1] + B[2]*S[2];

    double Q = B2 / D;
    double psi = BS*BS / (D*D*D);
    double e = tau/D;

    double n = (gamma_law-1.0)/gamma_law;

    //Newton-Raphson using the 2D scheme of Noble et al.
    //TODO: Rearrange to optimize for cold flows.

    //Initial guess: previous v2, eta
    double u[3] = {p[URR], p[UPP], p[UZZ]};
    double u2 = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    double w = sqrt(1.0 + u2);
    double v20 = u2 / (1.0 + u2);
    double eta0 = w*(1.0 + p[PPP]/(n*p[RHO]));

    //Run Newton-Raphson
    double v2, v21, eta, eta1;
    v21 = v20;
    eta1 = eta0;
    
    int i = 0;
    int clean = -1;

    if(DEBUG4)
        printf("      HLL Ps: s2=%.12lg e=%.12lg Q=%.12lg psi=%.12lg\n", 
                s2, e, Q, psi);

    while(1)
    {
        v2 = v21;
        eta = eta1;
        
        w = 1.0/sqrt(1.0-v2);
        double fa = (eta+Q)*(eta+Q)*v2 - psi*(2*eta+Q)/(eta*eta) - s2;
        double fb = eta - n*(1-v2)*(eta-w) + 0.5*Q*(1+v2) - 0.5*psi/(eta*eta)
                    - e;
        double dfadv2 = (eta+Q)*(eta+Q);
        double dfadet = 2*(eta+Q)*(v2 + psi/(eta*eta*eta));
        double dfbdv2 = n * (eta-0.5*w) + 0.5*Q;
        double dfbdet = 1 - n*(1-v2) + psi/(eta*eta*eta);
        double detj = dfadv2*dfbdet - dfbdv2*dfadet;

        v21  = v2  - ( dfbdet*fa - dfadet*fb)/detj;
        eta1 = eta - (-dfbdv2*fa + dfadv2*fb)/detj;

        i++;

        if(v21 > 1.0)
            v21 = 0.5*(v2 + 1.0);
        if(v21 < 0.0)
            v21 = 0.5*v2;
        if(eta1 < 1.0)
            eta1 = 0.5*(eta+1.0);

        if(DEBUG4)
        {
            printf("      %d: v2=%.12lg v21=%.12lg eta=%.12lg eta1=%.12lg\n",
                    i,  v2, v21, eta, eta1);
            printf("         fa=%.12lg fb=%.12lg\n", fa, fb);
        }

        double err = (eta1-eta)/eta;

        if(fabs(err) < prec && clean < 0)
            clean = Nextra+1;
        if(clean >= 0)
            clean--;
        if(clean == 0 || i == max_iter)
            break;
    }

    v2 = v21;
    eta = eta1;

    //Prim recovery
    w = 1.0/sqrt(1-v2);
    double h = eta/w;

    double rho = D / w;
    double Pp = n * rho * (h-1.0);
    
    h = 1.0 + gamma_law/(gamma_law-1.0) * Pp/rho;

    u[0] = (S[0] + BS*B[0]/(D*h*w)) / (D*h + B2/w);
    u[1] = (S[1] + BS*B[1]/(D*h*w)) / (D*h + B2/w);
    u[2] = (S[2] + BS*B[2]/(D*h*w)) / (D*h + B2/w);

    double uB = u[0]*B[0] + u[1]*B[1] + u[2]*B[2];
    double b2 = (B2 + uB*uB) / (w*w);

    p[RHO] = rho;
    p[URR] = u[0];
    p[UPP] = u[1];
    p[UZZ] = u[2];
    p[PPP] = Pp;

    return Pp + 0.5*b2;
}

