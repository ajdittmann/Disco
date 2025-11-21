#include "../paul.h"

static double q_planet = 1.0;
static double Npl = 3;

void rhs(int N, double eps, double *sim, double *update);
void rk8(int N, double dt, double *input, double *output);

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = Npl;
   q_planet = theDomain->theParList.Mass_Ratio;

}

int planet_motion_analytic( void ){
   return(1);
}

void initializePlanets( struct planet * thePlanets ){

   double a  = 1.0;

   double q = 1.0; //q_planet;
   double mu = q/(1.+q);

   double om = pow( a , -1.5 );

   thePlanets[0].M     = (1-mu);
   thePlanets[0].vr    = 0.0;
   thePlanets[0].omega = om;
   thePlanets[0].vz    = 0.0;
   thePlanets[0].r     = a*mu;
   thePlanets[0].phi   = M_PI;
   thePlanets[0].z     = 0.0;
   thePlanets[0].eps   = 0.05; //eps
   thePlanets[0].type  = PLPOINTMASS;

   thePlanets[1].M     = mu;
   thePlanets[1].vr    = 0.0;
   thePlanets[1].omega = om;
   thePlanets[1].vz    = 0.0;
   thePlanets[1].r     = a*(1.-mu);
   thePlanets[1].phi   = 0.0;
   thePlanets[1].z     = 0.0;
   thePlanets[1].eps   = 0.05; //eps;
   thePlanets[1].type  = PLPOINTMASS;


   thePlanets[2].M     = 0.0; //q_planet;
   thePlanets[2].vr    = 0.0;
   thePlanets[2].omega = pow(3.5, -1.5);
   thePlanets[2].vz    = 0.0;
   thePlanets[2].r     = 3.5;
   thePlanets[2].phi   = 0.0;
   thePlanets[2].z     = 0.0;
   thePlanets[2].eps   = 0.1;
   thePlanets[2].type  = PLPOINTMASS;
}

void movePlanets( struct planet * thePlanets , double t , double dt ){

  double *planetData = calloc(Npl*7, sizeof(double));
  double *planetDataOut = calloc(Npl*7, sizeof(double));

  double cp, sp;
  for (int k=0; k<Npl; k++){
    planetData[k*7] = thePlanets[k].M;
    cp = cos(thePlanets[k].phi); sp = sin(thePlanets[k].phi);
    planetData[k*7 + 1] = thePlanets[k].r*cp;
    planetData[k*7 + 2] = thePlanets[k].r*sp;
    planetData[k*7 + 3] = 0.0;
    planetData[k*7 + 4] = thePlanets[k].vr*cp - thePlanets[k].omega*thePlanets[k].r*sp;
    planetData[k*7 + 5] = thePlanets[k].vr*sp + thePlanets[k].omega*thePlanets[k].r*cp;
    planetData[k*7 + 6] = 0.0;
  }

  rk8(Npl, dt, planetData, planetDataOut);

  for (int k=0; k<Npl; k++){
    thePlanets[k].r = sqrt( SQR(planetDataOut[k*7+1]) + SQR(planetDataOut[k*7+2]));
    thePlanets[k].phi = atan2(planetDataOut[k*7+2], planetDataOut[k*7+1]);
    cp = cos(thePlanets[k].phi); sp = sin(thePlanets[k].phi);
    thePlanets[k].vr = planetDataOut[k*7+4]*cp + planetDataOut[k*7+5]*sp;
    thePlanets[k].omega = (-planetDataOut[k*7+4]*sp + planetDataOut[k*7+5]*cp)/thePlanets[k].r;
  }

  free(planetData);
  free(planetDataOut);
}



#define a2_0 0.25

#define a3_0  .8740084650491523205268632759487741197705e-1
#define a3_1  .2548760493865432175308795062034568513581e-1

#define a4_0  .4233316929133858267716535433070866141732e-1
#define a4_2  .1269995078740157480314960629921259842520

#define a5_0  .4260950588874226149488144523757227409094
#define a5_2 -1.598795284659152326542773323065718111709
#define a5_3  1.596700225771729711593958870689995370799

#define a6_0  .5071933729671392951509061813851363923933e-1
#define a6_3  .2543337726460040758275471440887777803137
#define a6_4  .2039468900572819946573622377727085804470

#define a7_0 -.2900037471752311097038837928542589612409
#define a7_3  1.344187391026078988943868110941433700318
#define a7_4 -2.864777943361442730961110382703656282947
#define a7_5  2.677594299510594851721126064616481543870

#define a8_0  .9853501133799354646974040298072701428476e-1
#define a8_4  .2219268063075138484202403649819738790358
#define a8_5 -.1814062291180699431269033828807395245747
#define a8_6  .1094441147256254823692261491803863125415e-1

#define a9_0  .3871105254573114467944461816516637340565
#define a9_3 -1.442445497485527757125674555307792776717
#define a9_4  2.905398189069950931769134644923384844174
#define a9_5 -1.853771069630105929084333267581197802518
#define a9_6  .1400364809872815426949732510977124147922
#define a9_7  .5727394081149581657574677462444770648875

#define a10_0 -.1612440344443930810063001619791348059544
#define a10_3 -.1733960295735898408357840447396256789490
#define a10_4 -1.301289281406514740601681274517249252974
#define a10_5  1.137950375173861730855879213143100347212
#define a10_6 -.3174764966396688010692352113804302469898e-1
#define a10_7  .9335129382493366643981106448605688485659
#define a10_8 -.8378631833473385270330085562961643320150e-1

#define a11_0 -.1919944488158953328151080465148357607314e-1
#define a11_3  .2733085726526428490794232625401612427562
#define a11_4 -.6753497320694437291969161121094238085624
#define a11_5  .3415184981384601607173848997472838271198
#define a11_6 -.6795006480337577247892051619852462939191e-1
#define a11_7  .9659175224762387888426558649121637650975e-1
#define a11_8  .1325308251118210118072103846654538995123
#define a11_9  .3685495936038611344690632995153166681295

#define a12_0  .6091877403645289867688841211158881778458
#define a12_3 -2.272569085898001676899980093141308839972
#define a12_4  4.757898342694029006815525588191478549755
#define a12_5 -5.516106706692758482429468966784424824484
#define a12_6  .2900596369680119270909581856594617437818
#define a12_7  .5691423963359036822910985845480184914563
#define a12_8  .7926795760332167027133991620589332757995
#define a12_9  .1547372045328882289412619077184989823205
#define a12_10  1.614970895662181624708321510633454443497

#define b0  .4472956466669571420301584042904938246647e-1
#define b5  .1569103352770819981336869801072664540918
#define b6  .1846097340815163774070245187352627789204
#define b7  .2251638060208699104247941940035072197092
#define b8  .1479461565197023468700517988544914175374
#define b9  .7605554244495582526979836191033649101273e-1
#define b10  .1227729023501861961082434631592143738854
#define b11  .4181195863899163158338484280087188237679e-1

void rk8(int N, double dt, double *input, double *output){
/*  Taken from:
#     J.H. Verner, SIAM NA 1978, 772-790,
#     Explicit Runge--Kutta methods with estimates of the
#     Local Truncation Error
*/
  int i;
  double *xn = calloc(N*7, sizeof(double));
  for (i=0; i<N*7; i++) xn[i] = input[i];

  double *k1 = calloc(N*7, sizeof(double));
  double *k2 = calloc(N*7, sizeof(double));
  double *k3 = calloc(N*7, sizeof(double));
  double *k4 = calloc(N*7, sizeof(double));
  double *k5 = calloc(N*7, sizeof(double));
  double *k6 = calloc(N*7, sizeof(double));
  double *k7 = calloc(N*7, sizeof(double));
  double *k8 = calloc(N*7, sizeof(double));
  double *k9 = calloc(N*7, sizeof(double));
  double *k10 = calloc(N*7, sizeof(double));
  double *k11 = calloc(N*7, sizeof(double));
  double *k12 = calloc(N*7, sizeof(double));

  rhs(N, 0.0, xn, k1);
  for (i=0; i<N*7; i++) xn[i] = input[i] + dt*a2_0*k1[i];
  rhs(N, 0.0, xn, k2);
  for (i=0; i<N*7; i++) xn[i] = input[i] + dt*(a3_0*k1[i] + a3_1*k2[i]);
  rhs(N, 0.0, xn, k3);
  for (i=0; i<N*7; i++) xn[i] = input[i] + dt*(a4_0*k1[i] + a4_2*k3[i]);
  rhs(N, 0.0, xn, k4);
  for (i=0; i<N*7; i++) xn[i] = input[i] + dt*(a5_0*k1[i] + a5_2*k3[i] + a5_3*k4[i]);
  rhs(N, 0.0, xn, k5);
  for (i=0; i<N*7; i++) xn[i] = input[i] + dt*(a6_0*k1[i] + a6_3*k4[i] + a6_4*k5[i]);
  rhs(N, 0.0, xn, k6);
  for (i=0; i<N*7; i++) xn[i] = input[i] + dt*(a7_0*k1[i] + a7_3*k4[i] + a7_4*k5[i] + a7_5*k6[i]);
  rhs(N, 0.0, xn, k7);
  for (i=0; i<N*7; i++) xn[i] = input[i] + dt*(a8_0*k1[i] + a8_4*k5[i] + a8_5*k6[i] + a8_6*k7[i]);
  rhs(N, 0.0, xn, k8);
  for (i=0; i<N*7; i++) xn[i] = input[i] + dt*(a9_0*k1[i] + a9_3*k4[i] + a9_4*k5[i] + a9_5*k6[i] + a9_6*k7[i] + a9_7*k8[i]);
  rhs(N, 0.0, xn, k9);
  for (i=0; i<N*7; i++) xn[i] = input[i] + dt*(a10_0*k1[i] + a10_3*k4[i] + a10_4*k5[i] + a10_5*k6[i] + a10_6*k7[i] + a10_7*k8[i] + a10_8*k9[i]);
  rhs(N, 0.0, xn, k10);
  for (i=0; i<N*7; i++) xn[i] = input[i] + dt*(a11_0*k1[i] + a11_3*k4[i] + a11_4*k5[i] + a11_5*k6[i] + a11_6*k7[i] + a11_7*k8[i] + a11_8*k9[i] + a11_9*k10[i]);
  rhs(N, 0.0, xn, k11);
  for (i=0; i<N*7; i++) xn[i] = input[i] + dt*(a12_0*k1[i] + a12_3*k4[i] + a12_4*k5[i] + a12_5*k6[i] + a12_6*k7[i] + a12_7*k8[i] + a12_8*k9[i] + a12_9*k10[i] + a12_10*k11[i]);
  rhs(N, 0.0, xn, k12);
  for (i=0; i<N*7; i++) output[i] = input[i] + dt*(b0*k1[i] + b5*k6[i] + b6*k7[i] + b7*k8[i] + b8*k9[i] + b9*k10[i] + b10*k11[i] + b11*k12[i]);

  free(xn);
  free(k1);
  free(k2);
  free(k3);
  free(k4);
  free(k5);
  free(k6);
  free(k7);
  free(k8);
  free(k9);
  free(k10);
  free(k11);
  free(k12);

}


void rhs(int N, double eps, double *sim, double *update){
  int i, j, k;
  double dx[3];
  double eps2 = eps*eps;
  double ri3, imag2, acci, rdr;
  for (i=0; i<N; i++){
    for (j=i+1; j<N; j++){
      rdr = 0.0;
      for (k=0; k<3; k++){
        dx[k] = sim[j*7 + 1 + k] - sim[i*7 + 1 + k];
        rdr += dx[k]*dx[k];
      }
      imag2 = 1.0/(rdr + eps2);
      ri3 = imag2*sqrt(imag2);
      for (k=0; k<3; k++){
        acci = ri3*dx[k];

        update[i*7+k+4] += acci*sim[j*7];
        update[j*7+k+4] -= acci*sim[i*7];
      }
    }
    for (k=0; k<3; k++){
      update[i*7+1+k] = sim[i*7+4+k];
    }
  }
}

