#include <cstdlib>
#include <cstring>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include "cir.hh"
#include <cmath>

/** Initializes the class for solving the CIR FP equation, 
 * using a variety of high-order and high-resolution methods.
 * \param[in] m_ the number of gridpoints to use.
 * \param[in] a_ the constant. 
 * \param[in] b_ the constant. 
 * \param[in] ss_ the constant. 
 * */
cir::cir(const int m_,const double ax_, const double bx_, 
        const double a_,const double b_,const double ss_) 
    : m(m_), dx(2./m), xsp(1/dx), ax(ax_), bx(bx_),
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

/** Computes the maximum timestep that can resolve the fluid advection, based
 * on the CFL condition.
 * \return The maximum timestep. */
double cir::advection_dt() {
    double adv_dt=0;
#pragma omp parallel for reduction(max:adv_dt)
    for(int j=0;j<m;j++) {
        // to modify rr depending on our domain
        double t, rr=j*dx;
        t=(a(b-rr)-ss)*xsp;
        if(t>adv_dt) adv_dt=t;
    }
    return adv_dt==0?std::numeric_limits<double>::max():1./adv_dt;
}

/** Chooses the timestep based on the limits from advection and viscosity.
 * \param[in] dt_pad the padding factor for the timestep for the physical
 *                   terms, which should be smaller than 1.
 * \param[in] adv_dt the maximum timestep to resolve the fluid advection.
 * \param[in] verbose whether to print out messages to the screen. */
void cir::choose_dt(double dt_pad,double adv_dt,bool verbose) {

    // Calculate the viscous timestep restriction
    double vis_dt=0.5*rho/(visc*(xxsp+yysp));
    int ca;

    // Choose the minimum of the two timestep restrictions
    if(adv_dt<vis_dt) {ca=0;dt_reg=adv_dt;}
    else {ca=1;dt_reg=vis_dt;}
    dt_reg*=dt_pad;

    // Print information if requested
    if(verbose) {
        const char mno[]="", myes[]=" <-- use this";
        printf("# Advection dt       : %g%s\n"
               "# Viscous dt         : %g%s\n"
               "# Padding factor     : %g\n"
               "# Minimum dt         : %g\n",
               adv_dt,ca==0?myes:mno,vis_dt,ca==1?myes:mno,dt_pad,dt_reg);
    }
}

/** Solves the transport equation, storing snapshots of the solution to a file.
 * \param[in] filename the name of the file to write to.
 * \param[in] snaps the number of snapshots to save (not including the initial
 *                  snapshot).
 * \param[in] duration the simulation duration.
 * \param[in] type the solve method to use. 0: FD, 1: FE, 2: DG (kiv) */
void cir::solve(const char* filename,int snaps,double duration,int type) {

    // Compute the timestep and number of iterations
    double adt;
    int iters=timestep_select(duration/snaps,adt);

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

double cir::true_density(double x_0, double x_t, double t_dur){
    double c = 2.*a/((1.-exp(-a*t_dur))*ss);
    double argu = x_t*2.*c;
    double dof = 4.*a*b/ss;
    double ncp = argu*exp(-a*t_dur);
    return pdf(non_central_chi_squared(dog,ncp),argu);
}