#include "GodonovSolver_t.h"
#include <cmath>

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

double GodonovSolver_t::u_star_initial_guess(double a_L, double a_R, double u_L, double u_R, double p_L, double p_R, double gamma_L, double gamma_R) {
    const double u_hat[] = {u_L + a_L * 2.0/(gamma_L - 1.0),
                            u_R - a_R * 2.0/(gamma_R - 1.0)};
    const double sigma = (p_L >= p_R) ? gamma_L : gamma_R;
    const double z = (gamma_L - 1.0)/(gamma_R - 1.0) * a_R/a_L * std::pow(p_L/p_R, 0.5 * (sigma - 1.0)/sigma);

    return (u_hat[0] * z + u_hat[1])/(1.0 + z);
}

void GodonovSolver_t::calculate_fluxes() {
    for (int i = 0; i < mesh_.n_faces_; ++i) {
        u_star_[i] = u_star_initial_guess(mesh_.a_[i], mesh_.a_[i+1], mesh_.u_[i], mesh_.u_[i+1], mesh_.p_[i], mesh_.p_[i+1], mesh_.gamma_[i], mesh_.gamma_[i+1]);
    }
}