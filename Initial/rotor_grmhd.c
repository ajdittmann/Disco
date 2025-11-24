#include "../disco.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

    double r   = x[0];
    double phi = x[1];

    double r0 = 0.1; 
    double Bx = 1.0;
    double omega_in = 9.995; 
    double omega_out = 0.0; 
    double Pp = 1.0;
    double rho_in = 10.0;
    double rho_out = 1.0;

    double rho = r > r0 ? rho_out : rho_in;
    double omega = r > r0 ? omega_out : omega_in;

    prim[RHO] = rho;
    prim[PPP] = Pp;
    prim[URR] = 0.0;
    prim[UPP] = omega*r*r / sqrt(1-r*r*omega*omega);
    prim[UZZ] = 0.0;

    prim[BRR] =  Bx*cos(phi);
    prim[BPP] = -Bx*sin(phi);
    prim[BZZ] = 0.0;

    if(NUM_N > 0)
    {
        int q;
        for(q = NUM_C; q < NUM_Q; q++)
            prim[q] = 0.0;
        
        if(r < r0)
            prim[NUM_C] = 1.0;
    }
}

