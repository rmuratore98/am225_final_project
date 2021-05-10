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
        /** The inverse grid spacing. */
        const double xsp;
        /** The lower bound on the domain. */
        const double ax;
        /** The upper bound on the domain. */
        const double bx;
        /** The constant, a. */
        const double a;
        /** The constant, b. */
        const double b;
        /** The constant, s^2. */
        const double ss;
        /** The regular timestep to be used. */
        double dt_reg;
        /** The primary grid for storing the solution. */
        double *p;
        /** The secondary grid for storing the solution. */
        double *q;
        /** The secondary grid for storing the advective solution. */
        double *t;
        /** The secondary grid for storing the diffusive solution. */
        double *d;
        cir(const int m_,const double ax_, const double bx_, const double a_,const double b_,const double ss_);
        ~cir();
        void initialize(double dt_pad,double max_vel=-1);
        void init_dirac();
        double advection_dt();
        void choose_dt(double dt_pad,double adv_dt,bool verbose=true);
        void solve(const char* filename,int snaps,double duration,int type);
        void fd(double dt);
        void fv(double dt);
        double true_density(double x_0, double x_t, double t_dur);
        double integral();
        /** Chooses a timestep size that is the largest value smaller than dt_reg,
        * such that a given interval length is a perfect multiple of this timestep.
        * \param[in] interval the interval length to consider.
        * \param[out] adt the timestep size.
        * \return The number of timesteps the fit into the interval. */
        inline int timestep_select(double interval, double &adt) {
            int l=static_cast<int>(interval/dt_reg)+1;
            adt=interval/l;
            return l;
        }
    private:
        inline double eno2(double p0,double p1,double p2,double p3);
        inline double min(double a,double b) {return a<b?a:b;}
        inline double max(double a,double b) {return a>b?a:b;}
        void print_line(FILE *fp,double x,double *zp,int snaps);
};

#endif