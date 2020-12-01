#include "GodonovSolver_t.h"

GodonovSolver_t::GodonovSolver_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double time, double discontinuity, int n_points, int problem_number) :
        Solver_t(rho_L, rho_R, u_L, u_R, p_L, p_R, time, discontinuity, n_points, problem_number){
    
}

GodonovSolver_t::~GodonovSolver_t() {}

void GodonovSolver_t::solve() {}