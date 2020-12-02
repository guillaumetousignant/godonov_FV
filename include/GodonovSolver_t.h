#ifndef GODONOVSOLVER_T_H
#define GODONOVSOLVER_T_H

#include "Solver_t.h"
#include "Mesh1D_t.h"
#include <vector>

template<typename FluxCalculator>
class GodonovSolver_t final : public Solver_t { 
    public: 
        GodonovSolver_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double time, double discontinuity, int n_points, int n_cells, int problem_number, double cfl);
        virtual ~GodonovSolver_t();

        Mesh1D_t mesh_;
        double cfl_;
        FluxCalculator flux_calculator_;

        virtual void solve();
        virtual void write_solution();

    private:
        double calculate_delta_t();
        void timestep(double delta_t);
};

#endif