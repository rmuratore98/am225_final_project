#ifndef CIR_HH
#define CIR_HH

#include <cstdio>
#include <cmath>

class cir {
    public:
        /** The grid size. */
        const int m;
        /** The grid spacing. */
        const double dx;
        /** The constant, a. */
        const double a;
        /** The constant, b. */
        const double b;
        /** The constant, s^2. */
        const double ss;
        /** The timestep. */
        double dt;
        /** The primary grid for storing the solution. */
        double *p;
        /** The secondary grid for storing the solution. */
        double *q;
        /** The secondary grid for storing the advective solution. */
        double *t;
        /** The secondary grid for storing the diffusive solution. */
        double *d;
        cir(int m_,double A_);
        ~cir();
        void init_dirac();
        void solve(const char* filename,int snaps,double duration,double safe_fac,int type);
        void fd(double dt);
        void fe(double dt);
        double integral();
    private:
        inline double eno2(double p0,double p1,double p2,double p3);
        void print_line(FILE *fp,double x,double *zp,int snaps);
};

#endif