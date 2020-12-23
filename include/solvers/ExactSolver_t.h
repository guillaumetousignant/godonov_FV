#ifndef FVM_EXACTSOLVER_T_H
#define FVM_EXACTSOLVER_T_H

#include "entities/Solver_t.h"
#include <vector>
#include <string>

namespace FVM { namespace Solvers {
    class ExactSolver_t final : public FVM::Entities::Solver_t { 
        public: 
            ExactSolver_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double x_L, double x_R, double time, double discontinuity, int n_points, int problem_number);
            virtual ~ExactSolver_t();

            double a_star_[2];
            double u_star_;
            double p_star_[2];
            double p_star_prime_[2];
            double C_[2];

            void left_shock();
            void right_shock();
            void left_rarefaction();
            void right_rarefaction();
            virtual void solve() override;
            void calculate_a_star();
            virtual void write_solution(std::string suffix = "") override;
    };
}}
#endif