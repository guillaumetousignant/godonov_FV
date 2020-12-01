#ifndef GODONOVSOLVER_T_H
#define GODONOVSOLVER_T_H

#include "Solver_t.h"
#include "Mesh1D_t.h"
#include <vector>

class GodonovSolver_t final : public Solver_t { 
    public: 
        GodonovSolver_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double time, double discontinuity, int n_cells, int problem_number);
        virtual ~GodonovSolver_t();

        Mesh1D_t mesh_;

        virtual void solve();
};

#endif