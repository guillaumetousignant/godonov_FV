#include "Solver_t.h"
#include <cmath>

Solver_t::Solver_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double time, double discontinuity, int n_cells, int problem_number) :
        gamma_{1.4, 1.4},
        R_(8.31446261815324),
        a_{std::sqrt(gamma_[0] * p_L/rho_L), std::sqrt(gamma_[1] * p_R/rho_R)},
        u_{u_L, u_R},
        p_{p_L, p_R},
        time_(time),
        discontinuity_(discontinuity),
        n_cells_(n_cells),
        problem_number_(problem_number),
        x_{0.0, 10.0} {


}

Solver_t::~Solver_t() {}