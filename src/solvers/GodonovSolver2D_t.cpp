#include "solvers/GodonovSolver2D_t.h"
#include "fluxes/ExactRiemannFlux_t.h"
#include "fluxes/RoeFlux_t.h"
#include "fluxes/RoeEntropyFlux_t.h"
#include "fluxes/HLLEFlux_t.h"
#include <cmath>
#include <algorithm>

using FVM::Solvers::GodonovSolver2D_t;
using FVM::Fluxes::ExactRiemannFlux_t;
using FVM::Fluxes::RoeFlux_t;
using FVM::Fluxes::RoeEntropyFlux_t;
using FVM::Fluxes::HLLEFlux_t;

template class GodonovSolver2D_t<ExactRiemannFlux_t>; // Like, I understand why I need this, but man is it crap.
template class GodonovSolver2D_t<RoeFlux_t>;
template class GodonovSolver2D_t<RoeEntropyFlux_t>;
template class GodonovSolver2D_t<HLLEFlux_t>;

template<typename FluxCalculator>
GodonovSolver2D_t<FluxCalculator>::GodonovSolver2D_t() :
        flux_calculator_() {}

template<typename FluxCalculator>
GodonovSolver2D_t<FluxCalculator>::~GodonovSolver2D_t() {}

template<typename FluxCalculator>
void GodonovSolver2D_t<FluxCalculator>::solve(double time, double cfl, FVM::Entities::Mesh2D_t& mesh) {
    double time = 0.0;

    while (time < end_time_) {
        double delta_t = calculate_delta_t();
        if (time + delta_t > end_time_) {
            delta_t = end_time_ - time;
        }

        mesh.boundary_conditions();
        flux_calculator_.calculate_fluxes(delta_t, mesh);
        timestep(delta_t, mesh_.delta_x_, mesh_.gamma_, mesh_.u_, mesh_.a_, mesh_.p_, mesh_.F_1_, mesh_.F_2_, mesh_.F_3_);
        time += delta_t;
    }
}

template<typename FluxCalculator>
void GodonovSolver2D_t<FluxCalculator>::write_solution(int problem_number, std::string suffix /* = "" */) {
    std::vector<double> x(n_points_ * mesh_.n_cells_);
    std::vector<double> rho(n_points_ * mesh_.n_cells_);
    std::vector<double> u(n_points_ * mesh_.n_cells_);
    std::vector<double> p(n_points_ * mesh_.n_cells_);
    std::vector<double> mach(n_points_ * mesh_.n_cells_);
    std::vector<double> T(n_points_ * mesh_.n_cells_);

    for (int i = 1; i <= mesh_.n_cells_; ++i) { // Don't take the ghost cells
        const int offset = (i - 1) * n_points_;
        for (int j = 0; j < n_points_; ++j) {
            // Watch out for integer division when lerping in c++
            x[offset + j] = (n_points_ == 1) ? mesh_.x_[i] : mesh_.x_[i] - mesh_.delta_x_/2 + j/(n_points_ - 1.0) * mesh_.delta_x_;
            rho[offset + j] = mesh_.gamma_[i] * mesh_.p_[i]/std::pow(mesh_.a_[i], 2);
            u[offset + j] = mesh_.u_[i];
            p[offset + j] = mesh_.p_[i];
            mach[offset + j] = u[offset + j] / std::sqrt(p[offset + j] / std::pow(rho[offset + j], mesh_.gamma_[i]));
            T[offset + j] = p[offset + j]/(rho[offset + j] * R_);
        }
    }

    write_file_data(n_points_ * mesh_.n_cells_, end_time_, rho, u, p, x, mach, T, problem_number_, suffix, mesh_.n_cells_);
}

template<typename FluxCalculator>
double GodonovSolver2D_t<FluxCalculator>::calculate_delta_t() {
    // Will have to do a reduce on a gpu?
    double max_u = 0.0;
    for (int i = 1; i <= mesh_.n_cells_; ++i) {
        max_u = std::max(max_u, std::abs(mesh_.u_[i]) + mesh_.a_[i]);
    }
    return cfl_ * mesh_.delta_x_/max_u;
}

template<typename FluxCalculator>
void GodonovSolver2D_t<FluxCalculator>::timestep(double delta_t, double delta_x, const std::vector<double> &gamma, std::vector<double> &u, std::vector<double> &a, std::vector<double> &p, const std::vector<double> F_1, const std::vector<double> F_2, const std::vector<double> F_3) {
    #pragma omp parallel for schedule(guided)
    for (int i = 1; i <= gamma.size() - 2; ++i) {
        double U_1 = gamma[i] * p[i] /std::pow(a[i], 2);
        double U_2 = gamma[i] * p[i] * u[i]/std::pow(a[i], 2);
        double U_3 = p[i]/(gamma[i] - 1) + gamma[i] * p[i] * std::pow(u[i], 2) * 0.5/std::pow(a[i], 2);

        U_1 += delta_t * (F_1[i-1] - F_1[i])/delta_x;
        U_2 += delta_t * (F_2[i-1] - F_2[i])/delta_x;
        U_3 += delta_t * (F_3[i-1] - F_3[i])/delta_x;

        u[i] = U_2/U_1;
        p[i] = (gamma[i] - 1) * (U_3 - U_2 * u[i] * 0.5);
        a[i] = std::sqrt(gamma[i] * p[i] /U_1);
    }
}