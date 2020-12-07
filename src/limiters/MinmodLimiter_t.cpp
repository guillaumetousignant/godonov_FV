#include "limiters/MinmodLimiter_t.h"
#include <algorithm>
#include <cmath>

FVM::Limiters::MinmodLimiter_t::MinmodLimiter_t() {}

FVM::Limiters::MinmodLimiter_t::~MinmodLimiter_t() {}

double FVM::Limiters::MinmodLimiter_t::phi(double a, double b) {
    return std::copysign(1.0, a) * std::max(0.0, std::min(std::abs(a), std::copysign(1.0, a) * b));
}
