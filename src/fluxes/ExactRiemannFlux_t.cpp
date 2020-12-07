#include "fluxes/ExactRiemannFlux_t.h"
#include "functions/RiemannProblem.h"
#include <cmath>

ExactRiemannFlux_t::ExactRiemannFlux_t(int n_faces) :
        a_star_L_(n_faces),
        a_star_R_(n_faces),
        u_star_(n_faces),
        p_star_L_(n_faces),
        p_star_R_(n_faces),
        p_star_prime_L_(n_faces),
        p_star_prime_R_(n_faces),
        C_L_(n_faces),
        C_R_(n_faces) {}

ExactRiemannFlux_t::~ExactRiemannFlux_t() {}

void ExactRiemannFlux_t::calculate_fluxes(double delta_t, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &F_1, std::vector<double> &F_2, std::vector<double> &F_3) {
    for (int i = 0; i < F_1.size(); ++i) {
        u_star_[i] = RiemannProblem::u_star_initial_guess(a[i], a[i+1], u[i], u[i+1], p[i], p[i+1], gamma[i], gamma[i+1]);
        C_L_[i] = RiemannProblem::calculate_C(a[i], p[i], gamma[i]);
        C_R_[i] = RiemannProblem::calculate_C(a[i+1], p[i+1], gamma[i+1]);
        RiemannProblem::solve_flux(a[i], a[i+1], u[i], u[i+1], p[i], p[i+1], gamma[i], gamma[i+1], C_L_[i], C_R_[i], u_star_[i], a_star_L_[i], a_star_R_[i], p_star_L_[i], p_star_R_[i], p_star_prime_L_[i], p_star_prime_R_[i]);
        RiemannProblem::calculate_a_star(a[i], a[i+1], u[i], u[i+1], p[i], p[i+1], gamma[i], gamma[i+1], u_star_[i], a_star_L_[i], a_star_R_[i], p_star_L_[i], p_star_R_[i]);
        double a_face, u_face, p_face, gamma_face;
        RiemannProblem::get_boundary_state(a_face, u_face, p_face, gamma_face, delta_t, a[i], a[i+1], u[i], u[i+1], p[i], p[i+1], gamma[i], gamma[i+1], u_star_[i], a_star_L_[i], a_star_R_[i], p_star_L_[i], p_star_R_[i]);

        F_1[i] = gamma_face * p_face * u_face /(std::pow(a_face, 2));
        F_2[i] = gamma_face * p_face * std::pow(u_face, 2)/std::pow(a_face, 2) + p_face;
        F_3[i] = u_face * (gamma_face * p_face /(gamma_face - 1) + gamma_face * p_face * std::pow(u_face, 2) * 0.5 /std::pow(a_face, 2));
    }
}