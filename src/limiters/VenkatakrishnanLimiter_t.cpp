#include "limiters/VenkatakrishnanLimiter_t.h"
#include "entities/Cell_t.h"
#include "entities/Face_t.h"
#include <algorithm>
#include <cmath>
#include <limits>

FVM::Limiters::VenkatakrishnanLimiter_t::VenkatakrishnanLimiter_t() {}

FVM::Limiters::VenkatakrishnanLimiter_t::~VenkatakrishnanLimiter_t() {}

void FVM::Limiters::VenkatakrishnanLimiter_t::calculate_derivatives(FVM::Entities::Mesh2D_t &mesh) {
    //#pragma omp parallel for schedule(guided)
    for (long long i = 0; i < mesh.n_cells_; ++i) {
        FVM::Entities::Cell_t& cell = mesh.cells_[i];

        cell.phi_a_ = 1;
        cell.phi_u_ = {1, 1};
        cell.phi_p_ = 1;
        cell.phi_gamma_ = 1;

        double a_max = -std::numeric_limits<double>::infinity();
        double a_min = std::numeric_limits<double>::infinity();
        FVM::Entities::Vec2f u_max = FVM::Entities::Vec2f{-std::numeric_limits<double>::infinity()};
        FVM::Entities::Vec2f u_min = FVM::Entities::Vec2f{std::numeric_limits<double>::infinity()};
        double p_max = -std::numeric_limits<double>::infinity();
        double p_min = std::numeric_limits<double>::infinity();
        double gamma_max = -std::numeric_limits<double>::infinity();
        double gamma_min = std::numeric_limits<double>::infinity();

        for (auto cell_index: cell.cells_) {
            const FVM::Entities::Cell_t& cell_k = mesh.cells_[cell_index];

            a_max = std::max(a_max, cell_k.a_);
            a_min = std::min(a_min, cell_k.a_);
            u_max = u_max.max(cell_k.u_);
            u_min = u_min.min(cell_k.u_);
            p_max = std::max(p_max, cell_k.p_);
            p_min = std::min(p_min, cell_k.p_);
            gamma_max = std::max(gamma_max, cell_k.gamma_);
            gamma_min = std::min(gamma_min, cell_k.gamma_);
        }

        for (auto face_index: cell.faces_) {
            const FVM::Entities::Face_t& face = mesh.faces_[face_index];

            const FVM::Entities::Vec2f delta = face.center_ - cell.center_;
            const double a_l = cell.a_ + delta.dot(cell.a_derivative_);
            const double ux_l = cell.u_.x() + delta.dot(cell.ux_derivative_);
            const double uy_l = cell.u_.y() + delta.dot(cell.uy_derivative_);
            const double p_l = cell.p_ + delta.dot(cell.p_derivative_);
            const double gamma_l = cell.gamma_ + delta.dot(cell.gamma_derivative_);

            double y_a;
            if (a_l > cell.a_) {
                y_a = (a_max - cell.a_)/(a_l - cell.a_);
            else if (a_l < cell.a_) {
                y_a = (a_max - cell.a_)/(a_l - cell.a_);
            }
            else {
                y_a = 2; // Will give phi = 1
            }
            const double phi_a_l = (std::pow(y_a, 2) + 2 * y_a)/(std::pow(y_a, 2) + y_a + 2);
            cell.phi_a_ = std::min(cell.phi_a_, phi_a_l);





        }
    }
}