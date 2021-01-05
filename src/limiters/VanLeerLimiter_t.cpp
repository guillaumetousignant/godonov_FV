#include "limiters/VanLeerLimiter_t.h"

FVM::Limiters::VanLeerLimiter_t::VanLeerLimiter_t() {}

FVM::Limiters::VanLeerLimiter_t::~VanLeerLimiter_t() {}

void FVM::Limiters::VanLeerLimiter_t::calculate_derivatives(const std::vector<double> &x, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &du_dx, std::vector<double> &da_dx, std::vector<double> &dp_dx) {
    constexpr double epsilon = 1.0e-16;

    #pragma omp parallel for schedule(static)
    for (int i = 1; i <= x.size() - 2; ++i) {
        const double delta_x_L = x[i] - x[i-1];
        const double delta_x_R = x[i+1] - x[i];

        const double a_u = (u[i] - u[i-1])/(delta_x_L);
        const double b_u = (u[i+1] - u[i])/(delta_x_R);
        const double a_a = (a[i] - a[i-1])/(delta_x_L);
        const double b_a = (a[i+1] - a[i])/(delta_x_R);
        const double a_p = (p[i] - p[i-1])/(delta_x_L);
        const double b_p = (p[i+1] - p[i])/(delta_x_R);

        du_dx[i] = (std::abs(a_u * b_u) + a_u * b_u)/(a_u + b_u + epsilon);
        da_dx[i] = (std::abs(a_a * b_a) + a_a * b_a)/(a_a + b_a + epsilon);
        dp_dx[i] = (std::abs(a_p * b_p) + a_p * b_p)/(a_p + b_p + epsilon);
    }
}