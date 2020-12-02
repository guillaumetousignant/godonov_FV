#include "RiemannProblem.h"
#include <cmath>
#include <limits>

double RiemannProblem::u_star_initial_guess(double a_L, double a_R, double u_L, double u_R, double p_L, double p_R, double gamma_L, double gamma_R) {
    const double u_hat_L = u_L + a_L * 2.0/(gamma_L - 1.0);
    const double u_hat_R = u_R - a_R * 2.0/(gamma_R - 1.0);
    const double sigma = (p_L >= p_R) ? gamma_L : gamma_R;
    const double z = (gamma_L - 1.0)/(gamma_R - 1.0) * a_R/a_L * std::pow(p_L/p_R, 0.5 * (sigma - 1.0)/sigma);

    return (u_hat_L * z + u_hat_R)/(1.0 + z);
}

double RiemannProblem::calculate_C(double a, double p, double gamma) {
    return gamma * p/a;
}

void RiemannProblem::calculate_a_star(double a_L, double a_R, double u_L, double u_R, double p_L, double p_R, double gamma_L, double gamma_R, double u_star, double a_star_L, double a_star_R, double p_star_L, double p_star_R) {
    if (u_star <= u_L) { // left shock
        a_star_L = a_L * std::sqrt(((gamma_L + 1.0) + (gamma_L - 1.0) * p_star_L/p_L) / ((gamma_L + 1.0) + (gamma_L - 1.0) * p_L/p_star_L));
    }
    else { // left rarefaction
        a_star_L = a_L - 0.5 * (gamma_L - 1.0) * (u_star - u_L); // Can't be 100% sure it was calculated I think
    }

    if (u_star >= u_R) { // right shock
        a_star_R = a_R * std::sqrt(((gamma_R + 1.0) + (gamma_R - 1.0) * p_star_R/p_R) / ((gamma_R + 1.0) + (gamma_R - 1.0) * p_R/p_star_R));
    }
    else { // right rarefaction
        a_star_R = a_R + 0.5 * (gamma_R - 1.0) * (u_star - u_R); // Can't be 100% sure it was calculated I think
    }
}

void RiemannProblem::left_shock(double a_L, double u_L, double p_L, double gamma_L, double C_L, double u_star, double &p_star_L, double &p_star_prime_L) {
    const double W = 0.25 * (gamma_L + 1.0) * (u_star - u_L) /a_L - std::sqrt(1.0 + std::pow(0.25 * (gamma_L + 1.0) * (u_star - u_L) /a_L, 2));

    p_star_L = p_L + C_L * (u_star - u_L) * W;
    p_star_prime_L = 2.0 * C_L * std::pow(W, 3) / (1.0 + std::pow(W, 2));
}

void RiemannProblem::right_shock(double a_R, double u_R, double p_R, double gamma_R, double C_R, double u_star, double &p_star_R, double &p_star_prime_R) {
    const double W = 0.25 * (gamma_R + 1.0) * (u_star - u_R) /a_R + std::sqrt(1.0 + std::pow(0.25 * (gamma_R + 1.0) * (u_star - u_R) /a_R, 2));

    p_star_R = p_R + C_R * (u_star - u_R) * W;
    p_star_prime_R = 2.0 * C_R * std::pow(W, 3) / (1.0 + std::pow(W, 2));
}

void RiemannProblem::left_rarefaction(double a_L, double u_L, double p_L, double gamma_L, double u_star, double &a_star_L, double &p_star_L, double &p_star_prime_L) {
    a_star_L = a_L - 0.5 * (gamma_L - 1.0) * (u_star - u_L);
    p_star_L = p_L * std::pow(a_star_L/a_L, 2 * gamma_L / (gamma_L - 1.0));
    p_star_prime_L = -gamma_L * p_star_L / a_star_L;
}

void RiemannProblem::right_rarefaction(double a_R, double u_R, double p_R, double gamma_R, double u_star, double &a_star_R, double &p_star_R, double &p_star_prime_R) {
    a_star_R = a_R + 0.5 * (gamma_R - 1.0) * (u_star - u_R);
    p_star_R = p_R * std::pow(a_star_R/a_R, 2 * gamma_R / (gamma_R - 1.0));
    p_star_prime_R = gamma_R * p_star_R / a_star_R;
}

void RiemannProblem::solve_flux(double a_L, double a_R, double u_L, double u_R, double p_L, double p_R, double gamma_L, double gamma_R, double C_L, double C_R, double &u_star, double a_star_L, double a_star_R, double &p_star_L, double &p_star_R, double &p_star_prime_L, double &p_star_prime_R) {
    constexpr double epsilon = 1.0e-6;
    double error = std::numeric_limits<double>::infinity();
    
    // This will always calculate at least once, since p_star_ is not calculated before
    do {
        (u_star <= u_L) ?  left_shock(a_L, u_L, p_L, gamma_L, C_L, u_star, p_star_L, p_star_prime_L) :  left_rarefaction(a_L, u_L, p_L, gamma_L, u_star, a_star_L, p_star_L, p_star_prime_L);
        (u_star >= u_R) ? right_shock(a_R, u_R, p_R, gamma_R, C_R, u_star, p_star_R, p_star_prime_R) : right_rarefaction(a_R, u_R, p_R, gamma_R, u_star, a_star_R, p_star_R,p_star_prime_R);

        u_star -= (p_star_L - p_star_R)/(p_star_prime_L - p_star_prime_R);
    }
    while (std::abs(1.0 - p_star_L/p_star_R) >= epsilon);
}