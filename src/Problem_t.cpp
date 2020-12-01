#include "Problem_t.h"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <sstream> 
#include <iomanip>
#include <limits>

namespace fs = std::filesystem;

Problem_t::Problem_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double time, double discontinuity, int n_points, int problem_number) :
        gamma_{1.4, 1.4},
        R_(8.31446261815324),
        a_{std::sqrt(gamma_[0] * p_L/rho_L), std::sqrt(gamma_[1] * p_R/rho_R)},
        u_{u_L, u_R},
        p_{p_L, p_R},
        time_(time),
        discontinuity_(discontinuity),
        n_points_(n_points),
        problem_number_(problem_number),
        x_span_{0.0, 10.0},
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

Problem_t::~Problem_t() {}

void Problem_t::left_shock() {
    const double W = 0.25 * (gamma_[0] + 1.0) * (u_star_ - u_[0]) /a_[0] - std::sqrt(1.0 + std::pow(0.25 * (gamma_[0] + 1.0) * (u_star_ - u_[0]) /a_[0], 2));

    p_star_[0] = p_[0] + C_[0] * (u_star_ - u_[0]) * W;
    p_star_prime_[0] = 2.0 * C_[0] * std::pow(W, 3) / (1.0 + std::pow(W, 2));
}

void Problem_t::right_shock() {
    const double W = 0.25 * (gamma_[1] + 1.0) * (u_star_ - u_[1]) /a_[1] + std::sqrt(1.0 + std::pow(0.25 * (gamma_[1] + 1.0) * (u_star_ - u_[1]) /a_[1], 2));

    p_star_[1] = p_[1] + C_[1] * (u_star_ - u_[1]) * W;
    p_star_prime_[1] = 2.0 * C_[1] * std::pow(W, 3) / (1.0 + std::pow(W, 2));
}

void Problem_t::left_rarefaction() {
    a_star_[0] = a_[0] - 0.5 * (gamma_[0] - 1.0) * (u_star_ - u_[0]);
    p_star_[0] = p_[0] * std::pow(a_star_[0]/a_[0], 2 * gamma_[0] / (gamma_[0] - 1.0));
    p_star_prime_[0] = -gamma_[0] * p_star_[0] / a_star_[0];
}

void Problem_t::right_rarefaction() {
    a_star_[1] = a_[1] + 0.5 * (gamma_[1] - 1.0) * (u_star_ - u_[1]);
    p_star_[1] = p_[1] * std::pow(a_star_[1]/a_[1], 2 * gamma_[1] / (gamma_[1] - 1.0));
    p_star_prime_[1] = gamma_[1] * p_star_[1] / a_star_[1];
}

void Problem_t::solve() {
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

void Problem_t::calculate_a_star() {
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

void Problem_t::write_file_data(int N_points, double time, const std::vector<double> &rho, const std::vector<double> &u, const std::vector<double> &p, const std::vector<double> &x, const std::vector<double> &mach, const std::vector<double> &T, int problem_number) {
    std::stringstream ss;
    std::ofstream file;

    fs::path save_dir = fs::current_path() / "data";
    fs::create_directory(save_dir);

    ss << "output_" << problem_number << ".dat";
    file.open(save_dir / ss.str());

    file << "TITLE = \"Problem " << problem_number << " at t= " << time << "\"" << std::endl;
    file << "VARIABLES = \"X\", \"U_x\", \"rho\", \"p\", \"mach\", \"T\"" << std::endl;
    file << "ZONE T= \"Zone     1\",  I= " << N_points << ",  J= 1,  DATAPACKING = POINT, SOLUTIONTIME = " << time << std::endl;

    for (int i = 0; i < N_points; ++i) {
        file << std::setw(12) << x[i] << " " << std::setw(12) << u[i] << " " << std::setw(12) << rho[i] << " " << std::setw(12) << p[i] << " " << std::setw(12) << mach[i] << " " << std::setw(12) << T[i] << std::endl;
    }

    file.close();
}

void Problem_t::write_solution() {
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

    const double wave_x[2][2] = {{discontinuity_ + time_ * wave_speed[0][0], discontinuity_ + time_ * wave_speed[0][1]}, 
                                {discontinuity_ + time_ * wave_speed[1][0], discontinuity_ + time_ * wave_speed[1][1]}};
    const double contact_x = discontinuity_ + time_ * u_star_;

    std::vector<double> x(n_points_);
    std::vector<double> rho(n_points_);
    std::vector<double> u(n_points_);
    std::vector<double> p(n_points_);
    std::vector<double> mach(n_points_);
    std::vector<double> T(n_points_);

    for (int i = 0; i < n_points_; ++i) {
        // Watch out for integer division when lerping in c++
        x[i] = x_span_[0] + i/(n_points_ - 1.0) * (x_span_[1] - x_span_[0]);

        if (x[i] < wave_x[0][0]) { // Left state
            rho[i] = gamma_[0] * p_[0]/std::pow(a_[0], 2);
            u[i] = u_[0];
            p[i] = p_[0];
            mach[i] = u[i] / std::sqrt(p[i] / std::pow(rho[i], gamma_[0]));
        }
        else if (x[i] < wave_x[0][1]) { // Left fan state
            const double v = (x[i] - discontinuity_)/time_;
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
            const double v = (x[i] - discontinuity_)/time_;
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

    write_file_data(n_points_, time_, rho, u, p, x, mach, T, problem_number_);
}