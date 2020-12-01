#include "GodonovSolver_t.h"
#include <cmath>

GodonovSolver_t::GodonovSolver_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double time, double discontinuity, int n_cells, int problem_number, double cfl) :
        Solver_t(rho_L, rho_R, u_L, u_R, p_L, p_R, time, discontinuity, n_cells, problem_number),
        mesh_(n_cells),
        cfl_(cfl) {

    mesh_.initial_conditions(a_[0], a_[1], u_[0], u_[1], p_[0], p_[1], x_[0], x_[1], discontinuity_);
}

GodonovSolver_t::~GodonovSolver_t() {}

void GodonovSolver_t::solve() {}

void GodonovSolver_t::calculate_fluxes() {

}