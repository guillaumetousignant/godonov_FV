#ifndef PROBLEM_T_H
#define PROBLEM_T_H

#include <vector>

class Problem_t { 
    public: 
        Problem_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double time, double discontinuity, int n_points, int problem_number);
        ~Problem_t();

        double gamma_[2];
        const double R_;
        double a_[2];
        double u_[2];
        double p_[2];
        double time_;
        double discontinuity_;
        int n_points_;
        int problem_number_;
        double x_span_[2];
        double a_star_[2];
        double u_star_;
        double p_star_[2];
        double p_star_prime_[2];
        double C_[2];

        void left_shock();
        void right_shock();
        void left_rarefaction();
        void right_rarefaction();
        void solve();
        void calculate_a_star();
        void write_file_data(int N_points, double time, const std::vector<double> &rho, const std::vector<double> &u, const std::vector<double> &p, const std::vector<double> &x, const std::vector<double> &mach, const std::vector<double> &T, int problem_number);
        void write_solution();
};

#endif