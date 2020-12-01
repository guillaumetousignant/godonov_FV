#include "Mesh1D_t.h"

Mesh1D_t::Mesh1D_t(int n_cells) :
        n_cells_(n_cells),
        rho_(new double[n_cells_ + 2]),  // Plus 2, because the ending cells are ghost cells
        u_(new double[n_cells_ + 2]),
        p_(new double[n_cells_ + 2]),
        x_(new double[n_cells_ + 2]),
        delta_x_((x_[1] - x_[0])/n_cells_),
        n_faces_(n_cells + 1) {}

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

void Mesh1D_t::initial_conditions(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double x_L, double x_R, double discontinuity) {
    for (int i = 0; i < n_cells_ + 2; ++i) { // No need for explicit treatment of the boundary conditions, nice
        x_[i] = x_L + (i - 0.5) * delta_x_;
        if (x_[i] < discontinuity) {
            rho_[i] = rho_L;
            u_[i] = u_L;
            p_[i] = p_L;
        }
        else {
            rho_[i] = rho_R;
            u_[i] = u_R;
            p_[i] = p_R;
        }
    }
}