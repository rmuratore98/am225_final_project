#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>

#include "cir.hh"

int main() {

    // Method: 0 = FD, 1 = FV;
    int type = 0;

    // Define other constants
    double a=2.,b=0.02,ss=0.15*0.15, x_0=0.8, bx=10., ax=0.;

    for(int i=10;i<40;i++){
        int m = int(10*pow(100,(1/30.)*i)+0.5);
        // Construct the simulation class, setting the number of gridpoints, the
        // periodicity, and physical constants
        cir cirfp(m,ax,bx,a,b,ss, x_0);
        //printf("dx = %g, xsp = %g\n", cirfp.dx, cirfp.xsp);
        //printf("%g\n", (bx-ax)/double(m));

        // Set the timestep based on multiplying the maximum allowable by a 
        // padding factor, initialize initial state
        cirfp.initialize(0.6,false);

        char buf[32];
        sprintf(buf,"cir_int_fd_%d.txt",1);

        // Run the simulation for a specified duration, outputting snapshots at
        // regular intervals
        double duration = 0.5;
        cirfp.solve(buf,10,duration,type);
        double q = cirfp.integral();
        printf("%d %g\n", m, q);
    
    }
    
}