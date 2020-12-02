#include "ExactRiemannFlux_t.h"
#include "RiemannProblem.h"
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

void ExactRiemannFlux_t::calculate_fluxes(Mesh1D_t &mesh, double delta_t) {
    for (int i = 0; i < mesh.n_faces_; ++i) {
        u_star_[i] = RiemannProblem::u_star_initial_guess(mesh.a_[i], mesh.a_[i+1], mesh.u_[i], mesh.u_[i+1], mesh.p_[i], mesh.p_[i+1], mesh.gamma_[i], mesh.gamma_[i+1]);
        C_L_[i] = RiemannProblem::calculate_C(mesh.a_[i], mesh.p_[i], mesh.gamma_[i]);
        C_R_[i] = RiemannProblem::calculate_C(mesh.a_[i+1], mesh.p_[i+1], mesh.gamma_[i+1]);
        RiemannProblem::solve_flux(mesh.a_[i], mesh.a_[i+1], mesh.u_[i], mesh.u_[i+1], mesh.p_[i], mesh.p_[i+1], mesh.gamma_[i], mesh.gamma_[i+1], C_L_[i], C_R_[i], u_star_[i], a_star_L_[i], a_star_R_[i], p_star_L_[i], p_star_R_[i], p_star_prime_L_[i], p_star_prime_R_[i]);
        RiemannProblem::calculate_a_star(mesh.a_[i], mesh.a_[i+1], mesh.u_[i], mesh.u_[i+1], mesh.p_[i], mesh.p_[i+1], mesh.gamma_[i], mesh.gamma_[i+1], u_star_[i], a_star_L_[i], a_star_R_[i], p_star_L_[i], p_star_R_[i]);
        double a_face, u_face, p_face, gamma_face;
        RiemannProblem::get_boundary_state(a_face, u_face, p_face, gamma_face, delta_t, mesh.a_[i], mesh.a_[i+1], mesh.u_[i], mesh.u_[i+1], mesh.p_[i], mesh.p_[i+1], mesh.gamma_[i], mesh.gamma_[i+1], u_star_[i], a_star_L_[i], a_star_R_[i], p_star_L_[i], p_star_R_[i]);

        mesh.F_1_[i] = gamma_face * p_face /(std::pow(a_face, 2));
        mesh.F_2_[i] = gamma_face * p_face * std::pow(u_face, 2)/std::pow(a_face, 2) + p_face;
        mesh.F_3_[i] = u_face * (gamma_face * p_face /(gamma_face - 1) + gamma_face * p_face * std::pow(u_face, 2) * 0.5 /std::pow(a_face, 2));
    }
}