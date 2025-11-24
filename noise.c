#include "disco.h"
#include "geometry.h"

static int type = 0;
static double absMag = 0.0;
static double relMag = 0.0;

void setNoiseParams(struct domain *theDomain)
{
    type = theDomain->theParList.noiseType;
    relMag = theDomain->theParList.noiseRel;
    absMag = theDomain->theParList.noiseAbs;
}

void addNoise(double *prim, double *x)
{
    if(type == 1)
    {
        double x1 = 2 * ((double)rand() / (double)(RAND_MAX)) - 1;
        double x2 = 2 * ((double)rand() / (double)(RAND_MAX)) - 1;
        double x3 = 2 * ((double)rand() / (double)(RAND_MAX)) - 1;
        double x4 = 2 * ((double)rand() / (double)(RAND_MAX)) - 1;

        prim[URR] = (1.0+x1*relMag)*prim[URR] + x2*absMag;
        prim[UPP] = (1.0+x3*relMag)*prim[UPP] + x4*absMag;
    }

    else if(type == 2)
    {
        double x1 = 2 * ((double)rand() / (double)(RAND_MAX)) - 1;
        double x2 = 2 * ((double)rand() / (double)(RAND_MAX)) - 1;

        prim[PPP] = (1.0+x1*relMag)*prim[PPP] + x2*absMag;
    }
    else if(type == 3)
    {
        double x1 = 2 * ((double)rand() / (double)(RAND_MAX)) - 1;
        double x2 = 2 * ((double)rand() / (double)(RAND_MAX)) - 1;
        double x3 = 2 * ((double)rand() / (double)(RAND_MAX)) - 1;
        double x4 = 2 * ((double)rand() / (double)(RAND_MAX)) - 1;

        double V[3] = {prim[URR], prim[UPP], prim[UZZ]};
        double Vc[3];
        get_vec_covariant(x, V, Vc); //Vc is now in orthonormal basis.

        prim[URR] = x1*relMag*Vc[1] + x2*absMag;
        prim[UZZ] = x3*relMag*Vc[1] + x4*absMag;
    }
}
