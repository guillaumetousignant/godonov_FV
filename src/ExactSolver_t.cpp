#include "ExactSolver_t.h"
#include <cmath>
#include <limits>

ExactSolver_t::ExactSolver_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double end_time, double discontinuity, int n_points, int problem_number) :
        Solver_t(rho_L, rho_R, u_L, u_R, p_L, p_R, end_time, discontinuity, n_points, problem_number),
        a_star_{0.0, 0.0},
        p_star_{0.0, 0.0},
        p_star_prime_{0.0, 0.0},
        C_{gamma_[0] * p_[0]/a_[0], gamma_[1] * p_[1]/a_[1]} {
    
    const double u_hat[] = {u_[0] + a_[0] * 2.0/(gamma_[0] - 1.0),
                            u_[1] - a_[1] * 2.0/(gamma_[1] - 1.0)};
    const double sigma = (p_[0] >= p_[1]) ? gamma_[0] : gamma_[1];
    const double z = (gamma_[0] - 1.0)/(gamma_[1] - 1.0) * a_[1]/a_[0] * std::pow(p_[0]/p_[1], 0.5 * (sigma - 1.0)/sigma);

    u_star_ = (u_hat[0] * z + u_hat[1])/(1.0 + z);
}

ExactSolver_t::~ExactSolver_t() {}

void ExactSolver_t::left_shock() {
    const double W = 0.25 * (gamma_[0] + 1.0) * (u_star_ - u_[0]) /a_[0] - std::sqrt(1.0 + std::pow(0.25 * (gamma_[0] + 1.0) * (u_star_ - u_[0]) /a_[0], 2));

    p_star_[0] = p_[0] + C_[0] * (u_star_ - u_[0]) * W;
    p_star_prime_[0] = 2.0 * C_[0] * std::pow(W, 3) / (1.0 + std::pow(W, 2));
}

void ExactSolver_t::right_shock() {
    const double W = 0.25 * (gamma_[1] + 1.0) * (u_star_ - u_[1]) /a_[1] + std::sqrt(1.0 + std::pow(0.25 * (gamma_[1] + 1.0) * (u_star_ - u_[1]) /a_[1], 2));

    p_star_[1] = p_[1] + C_[1] * (u_star_ - u_[1]) * W;
    p_star_prime_[1] = 2.0 * C_[1] * std::pow(W, 3) / (1.0 + std::pow(W, 2));
}

void ExactSolver_t::left_rarefaction() {
    a_star_[0] = a_[0] - 0.5 * (gamma_[0] - 1.0) * (u_star_ - u_[0]);
    p_star_[0] = p_[0] * std::pow(a_star_[0]/a_[0], 2 * gamma_[0] / (gamma_[0] - 1.0));
    p_star_prime_[0] = -gamma_[0] * p_star_[0] / a_star_[0];
}

void ExactSolver_t::right_rarefaction() {
    a_star_[1] = a_[1] + 0.5 * (gamma_[1] - 1.0) * (u_star_ - u_[1]);
    p_star_[1] = p_[1] * std::pow(a_star_[1]/a_[1], 2 * gamma_[1] / (gamma_[1] - 1.0));
    p_star_prime_[1] = gamma_[1] * p_star_[1] / a_star_[1];
}

void ExactSolver_t::solve() {
    constexpr double epsilon = 1.0e-6;
    double error = std::numeric_limits<double>::infinity();
    
    // This will always calculate at least once, since p_star_ is not calculated before
    do {
        (u_star_ <= u_[0]) ? left_shock() : left_rarefaction();
        (u_star_ >= u_[1]) ? right_shock() : right_rarefaction();

        u_star_ -= (p_star_[0] - p_star_[1])/(p_star_prime_[0] - p_star_prime_[1]);
    }
    while (std::abs(1.0 - p_star_[0]/p_star_[1]) >= epsilon);
}

void ExactSolver_t::calculate_a_star() {
    if (u_star_ <= u_[0]) { // left shock
        a_star_[0] = a_[0] * std::sqrt(((gamma_[0] + 1.0) + (gamma_[0] - 1.0) * p_star_[0]/p_[0]) / ((gamma_[0] + 1.0) + (gamma_[0] - 1.0) * p_[0]/p_star_[0]));
    }
    else { // left rarefaction
        a_star_[0] = a_[0] - 0.5 * (gamma_[0] - 1.0) * (u_star_ - u_[0]); // Can't be 100% sure it was calculated I think
    }

    if (u_star_ >= u_[1]) { // right shock
        a_star_[1] = a_[1] * std::sqrt(((gamma_[1] + 1.0) + (gamma_[1] - 1.0) * p_star_[1]/p_[1]) / ((gamma_[1] + 1.0) + (gamma_[1] - 1.0) * p_[1]/p_star_[1]));
    }
    else { // right rarefaction
        a_star_[1] = a_[1] + 0.5 * (gamma_[1] - 1.0) * (u_star_ - u_[1]); // Can't be 100% sure it was calculated I think
    }
}

void ExactSolver_t::write_solution(std::string suffix = "") {
    calculate_a_star(); // Can't be sure it was calculated

    double wave_speed[2][2]; // {{left_start, left_end}, {right_start, right_end}}

    if (u_star_ <= u_[0]) { // left shock
        const double mach = 0.25 * (gamma_[0] + 1.0) * (u_star_ - u_[0])/a_[0] - std::sqrt(1.0 + std::pow(0.25 * (gamma_[0] + 1.0) * (u_star_ - u_[0])/a_[0], 2));

        wave_speed[0][0] = wave_speed[0][1] = u_[0] + a_[0] * mach;
    }
    else { // left rarefaction
        wave_speed[0][0] = u_[0] - a_[0];
        wave_speed[0][1] = u_star_ - a_star_[0];
    }

    if (u_star_ >= u_[1]) { // right shock
        const double mach = 0.25 * (gamma_[1] + 1.0) * (u_star_ - u_[1])/a_[1] + std::sqrt(1.0 + std::pow(0.25 * (gamma_[1] + 1.0) * (u_star_ - u_[1])/a_[1], 2));

        wave_speed[1][0] = wave_speed[1][1] = u_[1] + a_[1] * mach;
    }
    else { // right rarefaction
        wave_speed[1][0] = u_star_ + a_star_[1];
        wave_speed[1][1] = u_[1] + a_[1];
    }

    const double wave_x[2][2] = {{discontinuity_ + end_time_ * wave_speed[0][0], discontinuity_ + end_time_ * wave_speed[0][1]}, 
                                {discontinuity_ + end_time_ * wave_speed[1][0], discontinuity_ + end_time_ * wave_speed[1][1]}};
    const double contact_x = discontinuity_ + end_time_ * u_star_;

    std::vector<double> x(n_points_);
    std::vector<double> rho(n_points_);
    std::vector<double> u(n_points_);
    std::vector<double> p(n_points_);
    std::vector<double> mach(n_points_);
    std::vector<double> T(n_points_);

    for (int i = 0; i < n_points_; ++i) {
        // Watch out for integer division when lerping in c++
        x[i] = x_[0] + i/(n_points_ - 1.0) * (x_[1] - x_[0]);

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

    write_file_data(n_points_, end_time_, rho, u, p, x, mach, T, problem_number_, suffix);
}