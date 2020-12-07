#include "limiters/SuperbeeLimiter_t.h"
#include <algorithm>
#include <cmath>


FVM::Limiters::SuperbeeLimiter_t::SuperbeeLimiter_t() {}

FVM::Limiters::SuperbeeLimiter_t::~SuperbeeLimiter_t() {}

double FVM::Limiters::SuperbeeLimiter_t::phi(double a, double b) {
    return std::copysign(1.0, a) * std::max(0.0, std::max(std::min(2 * std::abs(a), std::copysign(1.0, a) * b), std::min(std::abs(a), 2.0 * std::copysign(1.0, a) * std::abs(b))));
}