#include "limiters/VanLeerLimiter_t.h"

FVM::Limiters::VanLeerLimiter_t::VanLeerLimiter_t() {}

FVM::Limiters::VanLeerLimiter_t::~VanLeerLimiter_t() {}

double FVM::Limiters::VanLeerLimiter_t::phi(double a, double b) {
    constexpr double epsilon = 1.0e-16;
    return (std::abs(a * b) + a * b)/(a + b + epsilon);
}