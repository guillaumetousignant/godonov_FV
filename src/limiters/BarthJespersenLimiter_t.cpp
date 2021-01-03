#include "limiters/BarthJespersenLimiter_t.h"
#include "entities/Cell_t.h"
#include <algorithm>
#include <cmath>

FVM::Limiters::BarthJespersenLimiter_t::BarthJespersenLimiter_t() {}

FVM::Limiters::BarthJespersenLimiter_t::~BarthJespersenLimiter_t() {}

void FVM::Limiters::BarthJespersenLimiter_t::calculate_derivatives(FVM::Entities::Mesh2D_t &mesh) {
    //#pragma omp parallel for schedule(guided)
    for (long long i = 0; i < mesh.n_cells_; ++i) {
        FVM::Entities::Cell_t& cell = mesh.cells_[i];

    }
}

void FVM::Limiters::BarthJespersenLimiter_t::calculate_derivatives_hat(FVM::Entities::Mesh2D_t &mesh) {
    //#pragma omp parallel for schedule(guided)
    for (long long i = 0; i < mesh.n_cells_; ++i) {
        FVM::Entities::Cell_t& cell = mesh.cells_[i];

    }
}
