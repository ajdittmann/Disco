#include "../disco.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

    double rho, P, u;

    if(x[2] < 0.0)
    {
        rho = 1.0;
        P = 1.0;
        u = 0.0;
    }
    else
    {
        rho = 1.0;
        P = 0.1;
        u = 0.0;
    }


    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = 0.0;
    prim[UPP] = 0.0;
    prim[UZZ] = u;

    if(NUM_C >= 8)
    {
        prim[BRR] = 0.0;
        prim[BPP] = 0.0;
        prim[BZZ] = 1.0;
    }

    int q;
    for(q = NUM_C; q < NUM_Q; q++)
    {
        if(x[2] < 0.0)
            prim[q] = 1.0;
        else
            prim[q] = 0.0;
    }
}
