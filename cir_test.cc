#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>

#include "cir.hh"

const char fn[]="cir_fd.out";

int main() {

    // Create the output directory for storing the simulation frames
    mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Method: 0 = FD, 1 = FV;
    int type = 0;

    // Number of gridpoints
    const int m = 256;

    // Define other constants
    double a=2.,b=0.02,ss=0.025, x_0=0.8;

    // Construct the simulation class, setting the number of gridpoints, the
    // periodicity, and physical constants
    cir cirfp(m,0.,10.,a,b,ss, x_0);

    // Set the timestep based on multiplying the maximum allowable by a 
    // padding factor, initialize initial state
    cirfp.initialize(0.6);

    // Run the simulation for a specified duration, outputting snapshots at
    // regular intervals
    cirfp.solve(fn,10,.5,type);
    double q = cirfp.l2_loss(.5);
    printf("%i %g\n", m, q);
}