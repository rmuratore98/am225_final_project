#include <cstdlib>
#include <cstring>

#include "cir.hh"

/** Initializes the class for solving the CIR FP equation, 
 * using a variety of high-order and high-resolution methods.
 * \param[in] m_ the number of gridpoints to use.
 * \param[in] a_ the constant. 
 * \param[in] b_ the constant. 
 * \param[in] ss_ the constant. 
 * */
cir::cir(int m_,double a_,double b_,double ss_) : m(m_), dx(2./m), 
    a(a_), b(b_), ss(ss_), p(new double[m]), q(new double[m]) {}

/** The class destructor frees the dynamically allocated memory. */
cir::~cir() {
    delete [] q;
    delete [] p;
}

/** Initializes the solution to be a dirac delta function. */
void cir::init_dirac() {
    for(int i=0;i<m;i++) {
        double x=dx*(i+0.5);
        if(i==0) p[i]=1;
        else p[i]=0;
    }
}

/** Solves the transport equation, storing snapshots of the solution to a file.
 * \param[in] filename the name of the file to write to.
 * \param[in] snaps the number of snapshots to save (not including the initial
 *                  snapshot).
 * \param[in] duration the number of iterations to step the solution forward by
 *                     between snapshots.
 * \param[in] safe_fac a safety factor to apply to the CFL timestep restriction.
 * \param[in] type the solve method to use. 0: FD, 1: FE, 2: DG (kiv) */
void cir::solve(const char* filename,int snaps,double duration,double safe_fac,int type) {

    // Compute the timestep and number of iterations
    double interval=duration/snaps,dt=dx/A*safe_fac;
    int iters=static_cast<int>(interval/dt)+1;
    dt=interval/iters;

    // Allocate memory to store solution snapshots. Integrate the system and
    // store snapshots after fixed time intervals
    double *z=new double[m*(snaps+1)];
    memcpy(z,p,m*sizeof(double));
    for(int i=1;i<=snaps;i++) {

        // Perform the explict timesteps
        switch(type) {
            case 0: for(int k=0;k<iters;k++) fd(dt);break;
            case 1: for(int k=0;k<iters;k++) fe(dt);break;
        }

        // Store the snapshot
        memcpy(z+i*m,p,m*sizeof(double));
    }

    // Open the output file to store the snapshots
    FILE *fp=fopen(filename,"w");
    if(fp==NULL) {
        fputs("Can't open output file\n",stderr);
        exit(1);
    }

    // Print the snapshots, including periodic copies at either end to get a
    // full line over the interval from 0 to 1
    print_line(fp,-0.5*dx,z+(m-1),snaps);
    for(int j=0;j<m;j++) print_line(fp,(j+0.5)*dx,z+j,snaps);
    print_line(fp,1+0.5*dx,z,snaps);

    // Delete snapshots and close file
    fclose(fp);
    delete [] z;
}

/** Steps the simulation fields forward using FD.
 * \param[in] dt the time step to use. */
void cir::fd(dt){

    // Advective term
    double f=(a*(b-r)-ss)*2*dx;
    for(int j=0;j<m;j++) {
        // Compute required indices, taking into account periodicity (KIV)
        int jl=j==0?m-1+j:j-1,jll=j<=1?m-2+j:j-2,
            jr=j==m-1?1-m+j:j+1;

        t[j] = f*p[jr]-p[jl];

    }

    // Diffusive term
    double nu=dt/(dx*dx)*ss/2;
    for(int j=0;j<m;j++){
        // Compute indices on left and right, taking into account periodicity (KIV)
        int jl=j==0?m-1+j:j-1,
            jr=j==m-1?1-m+j:j+1;

        d[j]=p[j]+nu*(p[jl]-2*p[j]+p[jr]);
    }

    // Combining both terms
    for(int j=0;j<m;j++){
        q[j] = p[j] + dt*(-t[j]+d[j]+a*p[j]); 
    }

    // Swap pointers so that b becomes the primary array
    double *c=p;p=q;q=c;
}

/** Steps the simulation fields forward using FE.
 * \param[in] dt the time step to use. */
void cir::fe(dt){

    // Advective term
    double f=(a*(b-r)-ss)*2*dx;
    for(int j=0;j<m;j++) {

        // Compute required indices, taking into account periodicity (KIV)
        int jl=j==0?m-1+j:j-1,jll=j<=1?m-2+j:j-2,
            jr=j==m-1?1-m+j:j+1;

        // Perform update
        t[j]=p[j]-f*eno2(p[jr],p[j],p[jl],p[jll]);
    }

    // Diffusive term
    double nu=dt/(dx*dx)*ss/2;
    for(int j=0;j<m;j++){
        // Compute indices on left and right, taking into account periodicity (KIV)
        int jl=j==0?m-1+j:j-1,
            jr=j==m-1?1-m+j:j+1;

        d[j]=p[j]+nu*(p[jl]-2*p[j]+p[jr]);
    }

    // Combining both terms
    for(int j=0;j<m;j++){
        q[j] = p[j] + dt*(-t[j]+d[j]+a*p[j]); 
    }

    // Swap pointers so that b becomes the primary array
    double *c=p;p=q;q=c;
}

/** Calculates the ENO derivative using a sequence of values at four
 * gridpoints.
 * \param[in] (p0,p1,p2,p3) the sequence of values to use.
 * \return The computed derivative. */
inline double cir::eno2(double p0,double p1,double p2,double p3) {
    return fabs(p0-2*p1+p2)>fabs(p1-2*p2+p3)?3*p1-4*p2+p3:p0-p2;
}

/** Computes the integral of the solution over the domain.
 * \return The integral. */
double cir::integral() {
    double sum=0;
    for(double *ap=p;ap<p+m;ap++) sum+=*ap;
    return dx*sum;
}

/** Prints a line of stored snapshots to a file.
 * \param[in] fp a pointer to the file to write to.
 * \param[in] x the position in the domain corresponding to this line.
 * \param[in] zp a pointer to the first snapshot data point to print.
 * \param[in] snaps the number of snapshots (not including the starting
 *                  snapshot). */
void cir::print_line(FILE *fp,double x,double *zp,int snaps) {
    fprintf(fp,"%g",x);
    for(int i=0;i<=snaps;i++) fprintf(fp," %g",zp[i*m]);
    fputc('\n',fp);
}