#include "Mesh1D_t.h"

Mesh1D_t::Mesh1D_t(int n_cells) :
        n_cells_(n_cells),
        rho_(new double[n_cells_]), 
        u_(new double[n_cells_]),
        p_(new double[n_cells_]),
        x_(new double[n_cells_]) {}

Mesh1D_t::~Mesh1D_t() {
    if (rho_ != nullptr) {
        delete [] rho_;
    }

    if (u_ != nullptr) {
        delete [] u_;
    }

    if (p_ != nullptr) {
        delete [] p_;
    }

    if (x_ != nullptr) {
        delete [] x_;
    }
}