#ifndef FVM_GODONOVSOLVER_T_H
#define FVM_GODONOVSOLVER_T_H

#include "entities/Solver_t.h"
#include "entities/Mesh1D_t.h"
#include <vector>
#include <string>

template<typename FluxCalculator>
class GodonovSolver_t final : public Solver_t { 
    public: 
        GodonovSolver_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double x_L, double x_R, double time, double discontinuity, int n_points, int n_cells, int problem_number, double cfl);
        virtual ~GodonovSolver_t();

        Mesh1D_t mesh_;
        double cfl_;
        FluxCalculator flux_calculator_;

        virtual void solve();
        virtual void write_solution(std::string suffix = "");

    private:
        double calculate_delta_t();
        void timestep(double delta_t, double delta_x, const std::vector<double> &gamma, std::vector<double> &u, std::vector<double> &a, std::vector<double> &p, const std::vector<double> F_1, const std::vector<double> F_2, const std::vector<double> F_3);
};

#endif