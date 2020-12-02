#include "GodonovSolver_t.h"
#include "ExactRiemannFlux_t.h"
#include <cmath>
#include <algorithm>

template class GodonovSolver_t<ExactRiemannFlux_t>; // Like, I understand why I need this, but man is it crap.

template<typename FluxCalculator>
GodonovSolver_t<FluxCalculator>::GodonovSolver_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double x_L, double x_R, double time, double discontinuity, int n_points, int n_cells, int problem_number, double cfl) :
        Solver_t(rho_L, rho_R, u_L, u_R, p_L, p_R, x_L, x_R, time, discontinuity, n_points, problem_number),
        mesh_(n_cells, (x_[1] - x_[0])/n_cells),
        cfl_(cfl),
        flux_calculator_(n_cells + 1) {

    mesh_.initial_conditions(a_[0], a_[1], u_[0], u_[1], p_[0], p_[1], x_[0], x_[1], gamma_[0], gamma_[1], discontinuity_);
}

template<typename FluxCalculator>
GodonovSolver_t<FluxCalculator>::~GodonovSolver_t() {}

template<typename FluxCalculator>
void GodonovSolver_t<FluxCalculator>::solve() {
    double time = 0.0;

    while (time < end_time_) {
        double delta_t = calculate_delta_t();
        if (time + delta_t > end_time_) {
            delta_t = end_time_ - time;
        }

        flux_calculator_.calculate_fluxes(mesh_, delta_t);
        timestep(delta_t);
        time += delta_t;
    }
}

template<typename FluxCalculator>
void GodonovSolver_t<FluxCalculator>::write_solution(std::string suffix /* = "" */) {
    std::vector<double> x(n_points_ * mesh_.n_cells_);
    std::vector<double> rho(n_points_ * mesh_.n_cells_);
    std::vector<double> u(n_points_ * mesh_.n_cells_);
    std::vector<double> p(n_points_ * mesh_.n_cells_);
    std::vector<double> mach(n_points_ * mesh_.n_cells_);
    std::vector<double> T(n_points_ * mesh_.n_cells_);

    for (int i = 0; i < mesh_.n_cells_; ++i) {
        const int offset = i * n_points_;
        for (int j = 0; j < n_points_; ++j) {
            // Watch out for integer division when lerping in c++
            x[offset + j] = (n_points_ == 1) ? mesh_.x_[i] : mesh_.x_[i] - mesh_.delta_x_ + j/(n_points_ - 1.0) * mesh_.delta_x_;
            rho[offset + j] = mesh_.gamma_[i] * mesh_.p_[i]/std::pow(mesh_.a_[i], 2);
            u[offset + j] = mesh_.u_[i];
            p[offset + j] = mesh_.p_[i];
            mach[offset + j] = u[offset + j] / std::sqrt(p[offset + j] / std::pow(rho[offset + j], mesh_.gamma_[i]));
            T[offset + j] = p[offset + j]/(rho[offset + j] * R_);
        }
    }

    write_file_data(n_points_, end_time_, rho, u, p, x, mach, T, problem_number_, suffix);
}

template<typename FluxCalculator>
double GodonovSolver_t<FluxCalculator>::calculate_delta_t() {
    // Will have to do a reduce on a gpu?
    double max_u = 0.0;
    for (int i = 1; i <= mesh_.n_cells_; ++i) {
        max_u = std::max(max_u, std::abs(mesh_.u_[i]) + mesh_.a_[i]);
    }
    return cfl_ * mesh_.delta_x_/max_u;
}

template<typename FluxCalculator>
void GodonovSolver_t<FluxCalculator>::timestep(double delta_t) {
    for (int i = 1; i <= mesh_.n_cells_; ++i) {
        double U_1 = mesh_.gamma_[i] * mesh_.p_[i] /std::pow(mesh_.a_[i], 2);
        double U_2 = mesh_.gamma_[i] * mesh_.p_[i] * mesh_.u_[i]/std::pow(mesh_.a_[i], 2);
        double U_3 = mesh_.p_[i]/(mesh_.gamma_[i] - 1) + mesh_.gamma_[i] * mesh_.p_[i] * std::pow(mesh_.u_[i], 2) * 0.5/std::pow(mesh_.a_[i], 2);

        U_1 += delta_t * (mesh_.F_1_[i-1] - mesh_.F_1_[i])/mesh_.delta_x_;
        U_2 += delta_t * (mesh_.F_2_[i-1] - mesh_.F_2_[i])/mesh_.delta_x_;
        U_3 += delta_t * (mesh_.F_3_[i-1] - mesh_.F_3_[i])/mesh_.delta_x_;

        mesh_.u_[i] = U_2/U_1;
        mesh_.p_[i] = (mesh_.gamma_[i]) * (U_2 - U_1 * std::pow(mesh_.u_[i], 2) * 0.5);
        mesh_.a_[i] = std::sqrt(mesh_.gamma_[i] * mesh_.p_[i] /U_1);
    }
}