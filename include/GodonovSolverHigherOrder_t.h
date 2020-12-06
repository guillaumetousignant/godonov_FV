#ifndef GODONOVSOLVERHIGHERORDER_T_H
#define GODONOVSOLVERHIGHERORDER_T_H

#include "Solver_t.h"
#include "Mesh1D_t.h"
#include <vector>
#include <string>

template<typename FluxCalculator>
class GodonovSolverHigherOrder_t final : public Solver_t { 
    public: 
        GodonovSolverHigherOrder_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double x_L, double x_R, double time, double discontinuity, int n_points, int n_cells, int problem_number, double cfl);
        virtual ~GodonovSolverHigherOrder_t();

        Mesh1D_t mesh_;
        double cfl_;
        FluxCalculator flux_calculator_;
        std::vector<double> u_hat_;
        std::vector<double> a_hat_;
        std::vector<double> p_hat_;

        virtual void solve();
        virtual void write_solution(std::string suffix = "");

    private:
        double calculate_delta_t();
        void predictor(double delta_t);
        void timestep(double delta_t);
};

#endif