#include "entities/Mesh1D_t.h"

Mesh1D_t::Mesh1D_t(int n_cells, double delta_x) :
        n_cells_(n_cells),
        a_(n_cells_ + 2),  // Plus 2, because the ending cells are ghost cells
        u_(n_cells_ + 2),
        p_(n_cells_ + 2),
        x_(n_cells_ + 2),
        gamma_(n_cells_ + 2),
        delta_x_(delta_x),
        n_faces_(n_cells + 1),
        F_1_(n_faces_),
        F_2_(n_faces_),
        F_3_(n_faces_) {}

Mesh1D_t::~Mesh1D_t() {}

void Mesh1D_t::initial_conditions(double a_L, double a_R, double u_L, double u_R, double p_L, double p_R, double x_L, double x_R, double gamma_L, double gamma_R, double discontinuity) {
    for (int i = 0; i < n_cells_ + 2; ++i) { // No need for explicit treatment of the boundary conditions, nice
        x_[i] = x_L + (i - 0.5) * delta_x_;
        if (x_[i] < discontinuity) {
            a_[i] = a_L;
            u_[i] = u_L;
            p_[i] = p_L;
            gamma_[i] = gamma_L;
        }
        else {
            a_[i] = a_R;
            u_[i] = u_R;
            p_[i] = p_R;
            gamma_[i] = gamma_R;
        }
    }
}