#include "limiters/VanAlbadaLimiter_t.h"
#include <cmath>

FVM::Limiters::VanAlbadaLimiter_t::VanAlbadaLimiter_t() {}

FVM::Limiters::VanAlbadaLimiter_t::~VanAlbadaLimiter_t() {}

double FVM::Limiters::VanAlbadaLimiter_t::phi(double a, double b) {
    constexpr double epsilon = 1.0e-16;
    return (a * b * (a + b))/(std::pow(a, 2) + std::pow(b, 2) + epsilon);
}