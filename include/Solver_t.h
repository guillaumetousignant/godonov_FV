#ifndef SOLVER_T_H
#define SOLVER_T_H

class Solver_t { 
    public: 
        Solver_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double time, double discontinuity, int n_cells, int problem_number);
        virtual ~Solver_t();

        double gamma_[2];
        const double R_;
        double a_[2];
        double u_[2];
        double p_[2];
        double time_;
        double discontinuity_;
        int n_cells_;
        int problem_number_;
        double x_span_[2];

        virtual void solve() = 0;
};

#endif