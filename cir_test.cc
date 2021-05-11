#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>

#include "cir.hh"

//const char fn[]="cir_fd";
//const char file[]="cir_fd.out";

int main() {

    // Create the output directory for storing the simulation frames
    //mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Method: 0 = FD, 1 = FV;
    int type = 1;

    // Number of gridpoints
    const int m = 256;

    // Define other constants
    double a=2.,b=0.02,ss=0.15*0.15, x_0=0.8, bx=10., ax=0.;

    // Construct the simulation class, setting the number of gridpoints, the
    // periodicity, and physical constants
    cir cirfp(m,ax,bx,a,b,ss, x_0);
    //printf("dx = %g, xsp = %g\n", cirfp.dx, cirfp.xsp);
    //printf("%g\n", (bx-ax)/double(m));

    // Set the timestep based on multiplying the maximum allowable by a 
    // padding factor, initialize initial state
    cirfp.initialize(0.6,true);

    char buf[32];
    sprintf(buf,"cir_fv_%d.txt",1);

    // Run the simulation for a specified duration, outputting snapshots at
    // regular intervals
    double duration = 0.5;
    cirfp.solve(buf,10,duration,type,true);
    double q = cirfp.l2_loss(duration);
    printf("%i %g\n", m, q);
}