#include "GodonovSolver_t.h"
#include "RiemannProblem.h"

GodonovSolver_t::GodonovSolver_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double time, double discontinuity, int n_cells, int problem_number, double cfl) :
        Solver_t(rho_L, rho_R, u_L, u_R, p_L, p_R, time, discontinuity, n_cells, problem_number),
        mesh_(n_cells),
        cfl_(cfl),
        a_star_L_(new double[mesh_.n_faces_]),
        a_star_R_(new double[mesh_.n_faces_]),
        u_star_(new double[mesh_.n_faces_]),
        p_star_L_(new double[mesh_.n_faces_]),
        p_star_R_(new double[mesh_.n_faces_]),
        p_star_prime_L_(new double[mesh_.n_faces_]),
        p_star_prime_R_(new double[mesh_.n_faces_]),
        C_L_(new double[mesh_.n_faces_]),
        C_R_(new double[mesh_.n_faces_]) {

    mesh_.initial_conditions(a_[0], a_[1], u_[0], u_[1], p_[0], p_[1], x_[0], x_[1], gamma_[0], gamma_[1], discontinuity_);
}

GodonovSolver_t::~GodonovSolver_t() {
    if (a_star_L_ != nullptr) {
        delete [] a_star_L_;
    }

    if (a_star_R_ != nullptr) {
        delete [] a_star_R_;
    }

    if (u_star_ != nullptr) {
        delete [] u_star_;
    }

    if (p_star_L_ != nullptr) {
        delete [] p_star_L_;
    }

    if (p_star_R_ != nullptr) {
        delete [] p_star_R_;
    }

    if (p_star_prime_L_ != nullptr) {
        delete [] p_star_prime_L_;
    }

    if (p_star_prime_R_ != nullptr) {
        delete [] p_star_prime_R_;
    }

    if (C_L_ != nullptr) {
        delete [] C_L_;
    }

    if (C_R_ != nullptr) {
        delete [] C_R_;
    }
}

void GodonovSolver_t::solve() {}

void GodonovSolver_t::calculate_fluxes(double delta_t) {
    for (int i = 0; i < mesh_.n_faces_; ++i) {
        u_star_[i] = RiemannProblem::u_star_initial_guess(mesh_.a_[i], mesh_.a_[i+1], mesh_.u_[i], mesh_.u_[i+1], mesh_.p_[i], mesh_.p_[i+1], mesh_.gamma_[i], mesh_.gamma_[i+1]);
        C_L_[i] = RiemannProblem::calculate_C(mesh_.a_[i], mesh_.p_[i], mesh_.gamma_[i]);
        C_R_[i] = RiemannProblem::calculate_C(mesh_.a_[i+1], mesh_.p_[i+1], mesh_.gamma_[i+1]);
        RiemannProblem::solve_flux(mesh_.a_[i], mesh_.a_[i+1], mesh_.u_[i], mesh_.u_[i+1], mesh_.p_[i], mesh_.p_[i+1], mesh_.gamma_[i], mesh_.gamma_[i+1], C_L_[i], C_R_[i], u_star_[i], a_star_L_[i], a_star_R_[i], p_star_L_[i], p_star_R_[i], p_star_prime_L_[i], p_star_prime_R_[i]);
        RiemannProblem::calculate_a_star(mesh_.a_[i], mesh_.a_[i+1], mesh_.u_[i], mesh_.u_[i+1], mesh_.p_[i], mesh_.p_[i+1], mesh_.gamma_[i], mesh_.gamma_[i+1], u_star_[i], a_star_L_[i], a_star_R_[i], p_star_L_[i], p_star_R_[i]);
        double a_face, u_face, p_face;
        RiemannProblem::get_boundary_state(a_face, u_face, p_face);

    }
}