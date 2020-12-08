#ifndef FVM_GODONOVSOLVERHIGHERORDER_T_H
#define FVM_GODONOVSOLVERHIGHERORDER_T_H

#include "entities/Solver_t.h"
#include "entities/Mesh1D_t.h"
#include <vector>
#include <string>

namespace FVM { namespace Solvers {
    template<typename FluxCalculator, typename FluxLimiter>
    class GodonovSolverHigherOrder_t final : public FVM::Entities::Solver_t { 
        public: 
            GodonovSolverHigherOrder_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double x_L, double x_R, double time, double discontinuity, int n_points, int n_cells, int problem_number, double cfl);
            virtual ~GodonovSolverHigherOrder_t();

            FVM::Entities::Mesh1D_t mesh_;
            double cfl_;
            FluxCalculator flux_calculator_;
            FluxLimiter flux_limiter_;
            std::vector<double> u_hat_;
            std::vector<double> a_hat_;
            std::vector<double> p_hat_;
            std::vector<double> F_1_hat_;
            std::vector<double> F_2_hat_;
            std::vector<double> F_3_hat_;
            std::vector<double> du_dx_; // Limiter is included in this
            std::vector<double> da_dx_; // Limiter is included in this
            std::vector<double> dp_dx_hat_; // Limiter is included in this
            std::vector<double> du_dx_hat_; // Limiter is included in this
            std::vector<double> da_dx_hat_; // Limiter is included in this
            std::vector<double> dp_dx_hat_; // Limiter is included in this

            virtual void solve();
            virtual void write_solution(std::string suffix = "");

        private:
            double calculate_delta_t();
            void predictor(double delta_t, double delta_x, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &u_hat, std::vector<double> &a_hat, std::vector<double> &p_hat, const std::vector<double> F_1, const std::vector<double> F_2, const std::vector<double> F_3);
            void corrector(double delta_t, double delta_x, const std::vector<double> &gamma, std::vector<double> &u, std::vector<double> &a, std::vector<double> &p, const std::vector<double> F_1, const std::vector<double> F_2, const std::vector<double> F_3, const std::vector<double> F_1_hat, const std::vector<double> F_2_hat, const std::vector<double> F_3_hat);
    };
}}
#endif