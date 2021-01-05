#include "fluxes/ExactRiemannFlux_t.h"
#include "functions/RiemannProblem.h"
#include <cmath>

FVM::Fluxes::ExactRiemannFlux_t::ExactRiemannFlux_t() {}

FVM::Fluxes::ExactRiemannFlux_t::~ExactRiemannFlux_t() {}

void FVM::Fluxes::ExactRiemannFlux_t::calculate_fluxes(double delta_t, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &F_1, std::vector<double> &F_2, std::vector<double> &F_3) {
    #pragma omp parallel for schedule(static)
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
    #pragma omp parallel for schedule(static)
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
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < mesh.faces_.size(); ++i) {
        FVM::Entities::Face_t& face = mesh.faces_[i];
        const FVM::Entities::Cell_t& cell_L = mesh.cells_[face.cells_[0]];
        const FVM::Entities::Cell_t& cell_R = mesh.cells_[face.cells_[1]];

        const FVM::Entities::Vec2f u_prime_L(cell_L.u_.dot(face.normal_), cell_L.u_.dot(face.tangent_));
        const FVM::Entities::Vec2f u_prime_R(cell_R.u_.dot(face.normal_), cell_R.u_.dot(face.tangent_));

        double a_star_L, a_star_R;
        double p_star_L, p_star_R;

        double u_star = FVM::RiemannProblem::u_star_initial_guess(cell_L.a_, cell_R.a_, u_prime_L.x(), u_prime_R.x(), cell_L.p_, cell_R.p_, cell_L.gamma_, cell_R.gamma_);
        const double C_L = FVM::RiemannProblem::calculate_C(cell_L.a_, cell_L.p_, cell_L.gamma_);
        const double C_R = FVM::RiemannProblem::calculate_C(cell_R.a_, cell_R.p_, cell_R.gamma_);
        FVM::RiemannProblem::solve_flux(cell_L.a_, cell_R.a_, u_prime_L.x(), u_prime_R.x(), cell_L.p_, cell_R.p_, cell_L.gamma_, cell_R.gamma_, C_L, C_R, u_star, a_star_L, a_star_R, p_star_L, p_star_R);
        FVM::RiemannProblem::calculate_a_star(cell_L.a_, cell_R.a_, u_prime_L.x(), u_prime_R.x(), cell_L.p_, cell_R.p_, cell_L.gamma_, cell_R.gamma_, u_star, a_star_L, a_star_R, p_star_L, p_star_R);
        double a_face, u_face, p_face, gamma_face;
        FVM::RiemannProblem::get_boundary_state(a_face, u_face, p_face, gamma_face, delta_t, cell_L.a_, cell_R.a_, u_prime_L.x(), u_prime_R.x(), cell_L.p_, cell_R.p_, cell_L.gamma_, cell_R.gamma_, u_star, a_star_L, a_star_R, p_star_L, p_star_R);

        const FVM::Entities::Vec2f uv_face(u_face * face.normal_.x(), u_face * face.tangent_.x());
        face.F_1_ = u_face * gamma_face * p_face /(std::pow(a_face, 2)); // Those should probably be cached somewhere, they are computed twice.
        face.F_2_ = gamma_face * p_face * std::pow(uv_face.x(), 2)/std::pow(a_face, 2) + p_face;
        face.F_3_ = uv_face.x() * uv_face.y() * gamma_face * p_face/std::pow(a_face, 2);
        face.F_4_ = u_face * (gamma_face * p_face /(gamma_face - 1) + gamma_face * p_face * (std::pow(uv_face.x(), 2) + std::pow(uv_face.y(), 2)) * 0.5 /std::pow(a_face, 2));
    }
}

void FVM::Fluxes::ExactRiemannFlux_t::calculate_fluxes_higher_order(double delta_t, FVM::Entities::Mesh2D_t &mesh) {

}

void FVM::Fluxes::ExactRiemannFlux_t::calculate_fluxes_higher_order_hat(double delta_t, FVM::Entities::Mesh2D_t &mesh) {

}
