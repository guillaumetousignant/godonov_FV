#ifndef SOLVER_T_H
#define SOLVER_T_H

#include <vector>
#include <string>

class Solver_t { 
    public: 
        Solver_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double x_L, double x_R, double end_time, double discontinuity, int n_points, int problem_number);
        virtual ~Solver_t();

        double gamma_[2];
        const double R_;
        double a_[2];
        double u_[2];
        double p_[2];
        double end_time_;
        double discontinuity_;
        int n_points_; // Number of points per cell to plot
        int problem_number_;
        double x_[2];

        virtual void solve() = 0;
        virtual void write_solution(std::string suffix = "") = 0;
        void write_file_data(int N_points, double time, const std::vector<double> &rho, const std::vector<double> &u, const std::vector<double> &p, const std::vector<double> &x, const std::vector<double> &mach, const std::vector<double> &T, int problem_number, std::string suffix, int N);
};

#endif