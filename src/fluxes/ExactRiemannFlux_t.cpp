#include "fluxes/ExactRiemannFlux_t.h"
#include "functions/RiemannProblem.h"
#include <cmath>

FVM::Fluxes::ExactRiemannFlux_t::ExactRiemannFlux_t() {}

FVM::Fluxes::ExactRiemannFlux_t::~ExactRiemannFlux_t() {}

void FVM::Fluxes::ExactRiemannFlux_t::calculate_fluxes(double delta_t, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &F_1, std::vector<double> &F_2, std::vector<double> &F_3) {
    #pragma omp parallel for schedule(guided)
    for (int i = 0; i < F_1.size(); ++i) {
        double a_star_L, a_star_R;
        double p_star_L, p_star_R;

        double u_star = FVM::RiemannProblem::u_star_initial_guess(a[i], a[i+1], u[i], u[i+1], p[i], p[i+1], gamma[i], gamma[i+1]);
        const double C_L = FVM::RiemannProblem::calculate_C(a[i], p[i], gamma[i]);
        const double C_R = FVM::RiemannProblem::calculate_C(a[i+1], p[i+1], gamma[i+1]);
        FVM::RiemannProblem::solve_flux(a[i], a[i+1], u[i], u[i+1], p[i], p[i+1], gamma[i], gamma[i+1], C_L, C_R, u_star, a_star_L, a_star_R, p_star_L, p_star_R);
        FVM::RiemannProblem::calculate_a_star(a[i], a[i+1], u[i], u[i+1], p[i], p[i+1], gamma[i], gamma[i+1], u_star, a_star_L, a_star_R, p_star_L, p_star_R);
        double a_face, u_face, p_face, gamma_face;
        FVM::RiemannProblem::get_boundary_state(a_face, u_face, p_face, gamma_face, delta_t, a[i], a[i+1], u[i], u[i+1], p[i], p[i+1], gamma[i], gamma[i+1], u_star, a_star_L, a_star_R, p_star_L, p_star_R);

        F_1[i] = gamma_face * p_face * u_face /(std::pow(a_face, 2));
        F_2[i] = gamma_face * p_face * std::pow(u_face, 2)/std::pow(a_face, 2) + p_face;
        F_3[i] = u_face * (gamma_face * p_face /(gamma_face - 1) + gamma_face * p_face * std::pow(u_face, 2) * 0.5 /std::pow(a_face, 2));
    }
}

void FVM::Fluxes::ExactRiemannFlux_t::calculate_fluxes_higher_order(double delta_t, const std::vector<double> x, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &F_1, std::vector<double> &F_2, std::vector<double> &F_3, const std::vector<double> du_dx, const std::vector<double> da_dx, const std::vector<double> dp_dx) {
    #pragma omp parallel for schedule(guided)
    for (int i = 0; i < F_1.size(); ++i) {
        double a_star_L, a_star_R;
        double p_star_L, p_star_R;

        const double delta_x = (x[i+1] - x[i])/2; // CHECK this won't work with non-uniform meshes! Use (x - x_i)

        const double gamma_L = gamma[i];
        const double u_L = u[i] + du_dx[i] * delta_x;
        const double a_L = a[i] + da_dx[i] * delta_x;
        const double p_L = p[i] + dp_dx[i] * delta_x;
        const double gamma_R = gamma[i+1];
        const double u_R = u[i+1] - du_dx[i+1] * delta_x;
        const double a_R = a[i+1] - da_dx[i+1] * delta_x;
        const double p_R = p[i+1] - dp_dx[i+1] * delta_x;

        double u_star = FVM::RiemannProblem::u_star_initial_guess(a_L, a_R, u_L, u_R, p_L, p_R, gamma_L, gamma_R);
        const double C_L = FVM::RiemannProblem::calculate_C(a_L, p_L, gamma_L);
        const double C_R = FVM::RiemannProblem::calculate_C(a_R, p_R, gamma_R);
        FVM::RiemannProblem::solve_flux(a_L, a_R, u_L, u_R, p_L, p_R, gamma_L, gamma_R, C_L, C_R, u_star, a_star_L, a_star_R, p_star_L, p_star_R);
        FVM::RiemannProblem::calculate_a_star(a_L, a_R, u_L, u_R, p_L, p_R, gamma_L, gamma_R, u_star, a_star_L, a_star_R, p_star_L, p_star_R);
        double a_face, u_face, p_face, gamma_face;
        FVM::RiemannProblem::get_boundary_state(a_face, u_face, p_face, gamma_face, delta_t, a_L, a_R, u_L, u_R, p_L, p_R, gamma_L, gamma_R, u_star, a_star_L, a_star_R, p_star_L, p_star_R);

        F_1[i] = gamma_face * p_face * u_face /(std::pow(a_face, 2));
        F_2[i] = gamma_face * p_face * std::pow(u_face, 2)/std::pow(a_face, 2) + p_face;
        F_3[i] = u_face * (gamma_face * p_face /(gamma_face - 1) + gamma_face * p_face * std::pow(u_face, 2) * 0.5 /std::pow(a_face, 2));
    }
}

void FVM::Fluxes::ExactRiemannFlux_t::calculate_fluxes(double delta_t, FVM::Entities::Mesh2D_t &mesh) {
    
}
