#include "limiters/MinmodLimiter_t.h"
#include <algorithm>
#include <cmath>

FVM::Limiters::MinmodLimiter_t::MinmodLimiter_t() {}

FVM::Limiters::MinmodLimiter_t::~MinmodLimiter_t() {}

void FVM::Limiters::MinmodLimiter_t::calculate_derivatives(const std::vector<double> &x, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &du_dx, std::vector<double> &da_dx, std::vector<double> &dp_dx) {
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

        du_dx[i] = std::copysign(1.0, a_u) * std::max(0.0, std::min(std::abs(a_u), std::copysign(1.0, a_u) * b_u));
        da_dx[i] = std::copysign(1.0, a_a) * std::max(0.0, std::min(std::abs(a_a), std::copysign(1.0, a_a) * b_a));
        dp_dx[i] = std::copysign(1.0, a_p) * std::max(0.0, std::min(std::abs(a_p), std::copysign(1.0, a_p) * b_p));
    }
}
