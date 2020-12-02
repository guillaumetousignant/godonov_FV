#include "GodonovSolver_t.h"
#include <cmath>
#include <algorithm>

template<typename FluxCalculator>
GodonovSolver_t<FluxCalculator>::GodonovSolver_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double time, double discontinuity, int n_points, int n_cells, int problem_number, double cfl) :
        Solver_t(rho_L, rho_R, u_L, u_R, p_L, p_R, time, discontinuity, n_points, problem_number),
        mesh_(n_cells),
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
void GodonovSolver_t<FluxCalculator>::write_solution() {
    calculate_a_star(); // Can't be sure it was calculated

    std::vector<double> x(n_points_);
    std::vector<double> rho(n_points_);
    std::vector<double> u(n_points_);
    std::vector<double> p(n_points_);
    std::vector<double> mach(n_points_);
    std::vector<double> T(n_points_);

    for (int i = 0; i < mesh_.n_cells_; ++i) {
        // Watch out for integer division when lerping in c++
        x[i] = x_[0] + i/(mesh_.n_cells_ - 1.0) * (x_[1] - x_[0]);

        if (x[i] < wave_x[0][0]) { // Left state
            rho[i] = gamma_[0] * p_[0]/std::pow(a_[0], 2);
            u[i] = u_[0];
            p[i] = p_[0];
            mach[i] = u[i] / std::sqrt(p[i] / std::pow(rho[i], gamma_[0]));
        }
        else if (x[i] < wave_x[0][1]) { // Left fan state
            const double v = (x[i] - discontinuity_)/end_time_;
            const double a = (gamma_[0] - 1.0)/(gamma_[0] + 1.0) * (u_[0] - v) + 2.0/(gamma_[0] + 1.0) * a_[0];

            u[i] = v + a;
            p[i] = p_[0] * std::pow(a/a_[0], 2.0 * gamma_[0]/(gamma_[0] - 1.0));
            rho[i] = gamma_[0] * p[i] /std::pow(a, 2);
            mach[i] = u[i] / std::sqrt(p[i] / std::pow(rho[i], gamma_[0]));
        }
        else if (x[i] < contact_x) { // Left star state
            rho[i] = gamma_[0] * p_star_[0]/std::pow(a_star_[0], 2);
            u[i] = u_star_;
            p[i] = p_star_[0];
            mach[i] = u[i] / std::sqrt(p[i] / std::pow(rho[i], gamma_[0]));
        }
        else if (x[i] < wave_x[1][0]) { // Right star state
            rho[i] = gamma_[1] * p_star_[1]/std::pow(a_star_[1], 2);
            u[i] = u_star_;
            p[i] = p_star_[1];
            mach[i] = u[i] / std::sqrt(p[i] / std::pow(rho[i], gamma_[1]));
        }
        else if (x[i] < wave_x[1][1]) { // Right fan state
            const double v = (x[i] - discontinuity_)/end_time_;
            const double a = (gamma_[1] - 1.0)/(gamma_[1] + 1.0) * (v - u_[1]) + 2.0/(gamma_[1] + 1.0) * a_[1];

            u[i] = v - a;
            p[i] = p_[1] * std::pow(a/a_[1], 2.0 * gamma_[1]/(gamma_[1] - 1.0));
            rho[i] = gamma_[1] * p[i] /std::pow(a, 2);
            mach[i] = u[i] / std::sqrt(p[i] / std::pow(rho[i], gamma_[1]));
        }
        else { // Right state
            rho[i] = gamma_[1] * p_[1]/std::pow(a_[1], 2);
            u[i] = u_[1];
            p[i] = p_[1];
            mach[i] = u[i] / std::sqrt(p[i] / std::pow(rho[i], gamma_[1]));
        }

        T[i] = p[i]/(rho[i] * R_);
    }

    write_file_data(n_points_, end_time_, rho, u, p, x, mach, T, problem_number_);
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