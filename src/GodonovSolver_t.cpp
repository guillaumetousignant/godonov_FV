#include "GodonovSolver_t.h"

GodonovSolver_t::GodonovSolver_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double time, double discontinuity, int n_cells, int problem_number) :
        Solver_t(rho_L, rho_R, u_L, u_R, p_L, p_R, time, discontinuity, n_cells, problem_number),
        mesh_(n_cells) {

    for (int i = 0; i < n_cells_; ++i) {
        
    }
}

GodonovSolver_t::~GodonovSolver_t() {}

void GodonovSolver_t::solve() {}