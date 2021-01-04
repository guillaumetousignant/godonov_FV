#include "limiters/VenkatakrishnanLimiter_t.h"
#include "entities/Cell_t.h"
#include "entities/Face_t.h"
#include <algorithm>
#include <cmath>
#include <limits>

using FVM::Entities::Vec2f;

FVM::Limiters::VenkatakrishnanLimiter_t::VenkatakrishnanLimiter_t() {}

FVM::Limiters::VenkatakrishnanLimiter_t::~VenkatakrishnanLimiter_t() {}

void FVM::Limiters::VenkatakrishnanLimiter_t::calculate_derivatives(FVM::Entities::Mesh2D_t &mesh) {
    #pragma omp parallel for schedule(guided)
    for (long long i = 0; i < mesh.n_cells_; ++i) {
        FVM::Entities::Cell_t& cell = mesh.cells_[i];

        double phi_a = 1;
        Vec2f phi_u {1, 1};
        double phi_p = 1;
        double phi_gamma = 1;

        double a_max = -std::numeric_limits<double>::infinity();
        double a_min = std::numeric_limits<double>::infinity();
        Vec2f u_max {-std::numeric_limits<double>::infinity()};
        Vec2f u_min {std::numeric_limits<double>::infinity()};
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

            const Vec2f delta = face.center_ - cell.center_;
            const double a_l = cell.a_ + delta.dot(cell.a_derivative_);
            const Vec2f u_l = cell.u_ + Vec2f(delta.dot(cell.ux_derivative_),
                                              delta.dot(cell.uy_derivative_));
            const double p_l = cell.p_ + delta.dot(cell.p_derivative_);
            const double gamma_l = cell.gamma_ + delta.dot(cell.gamma_derivative_);

            double y_a;
            double y_ux;
            double y_uy;
            double y_p;
            double y_gamma;

            if (a_l > cell.a_) {
                y_a = (a_max - cell.a_)/(a_l - cell.a_);
            }
            else if (a_l < cell.a_) {
                y_a = (a_min - cell.a_)/(a_l - cell.a_);
            }
            else {
                y_a = 2; // Will give phi = 1
            }

            if (u_l.x() > cell.u_.x()) {
                y_ux = (u_max.x() - cell.u_.x())/(u_l.x() - cell.u_.x());
            }
            else if (u_l.x() < cell.u_.x()) {
                y_ux = (u_min.x() - cell.u_.x())/(u_l.x() - cell.u_.x());
            }
            else {
                y_ux = 2; // Will give phi = 1
            }

            if (u_l.y() > cell.u_.y()) {
                y_uy = (u_max.y() - cell.u_.y())/(u_l.y() - cell.u_.y());
            }
            else if (u_l.y() < cell.u_.y()) {
                y_uy = (u_min.y() - cell.u_.y())/(u_l.y() - cell.u_.y());
            }
            else {
                y_uy = 2; // Will give phi = 1
            }

            if (p_l > cell.p_) {
                y_p = (p_max - cell.p_)/(p_l - cell.p_);
            }
            else if (p_l < cell.p_) {
                y_p = (p_min - cell.p_)/(p_l - cell.p_);
            }
            else {
                y_p = 2; // Will give phi = 1
            }

            if (gamma_l > cell.gamma_) {
                y_gamma = (gamma_max - cell.gamma_)/(gamma_l - cell.gamma_);
            }
            else if (gamma_l < cell.gamma_) {
                y_gamma = (gamma_min - cell.gamma_)/(gamma_l - cell.gamma_);
            }
            else {
                y_gamma = 2; // Will give phi = 1
            }

            const double phi_a_l = (std::pow(y_a, 2) + 2 * y_a)/(std::pow(y_a, 2) + y_a + 2);
            const Vec2f phi_u_l {(std::pow(y_ux, 2) + 2 * y_ux)/(std::pow(y_ux, 2) + y_ux + 2),
                                 (std::pow(y_uy, 2) + 2 * y_uy)/(std::pow(y_uy, 2) + y_uy + 2)};
            const double phi_p_l = (std::pow(y_p, 2) + 2 * y_p)/(std::pow(y_p, 2) + y_p + 2);
            const double phi_gamma_l = (std::pow(y_gamma, 2) + 2 * y_gamma)/(std::pow(y_gamma, 2) + y_gamma + 2);

            phi_a = std::min(phi_a, phi_a_l);
            phi_u.min(phi_u_l);
            phi_p = std::min(phi_p, phi_p_l);
            phi_gamma = std::min(phi_gamma, phi_gamma_l);
        }

        cell.a_derivative_ *= phi_a;
        cell.ux_derivative_ *= phi_u.x();
        cell.uy_derivative_ *= phi_u.y();
        cell.p_derivative_ *= phi_p;
        cell.gamma_derivative_ *= phi_gamma;
    }
}

void FVM::Limiters::VenkatakrishnanLimiter_t::calculate_derivatives_hat(FVM::Entities::Mesh2D_t &mesh) {
    #pragma omp parallel for schedule(guided)
    for (long long i = 0; i < mesh.n_cells_; ++i) {
        FVM::Entities::Cell_t& cell = mesh.cells_[i];

        double phi_a = 1;
        Vec2f phi_u {1, 1};
        double phi_p = 1;
        double phi_gamma = 1;

        double a_max = -std::numeric_limits<double>::infinity();
        double a_min = std::numeric_limits<double>::infinity();
        Vec2f u_max {-std::numeric_limits<double>::infinity()};
        Vec2f u_min {std::numeric_limits<double>::infinity()};
        double p_max = -std::numeric_limits<double>::infinity();
        double p_min = std::numeric_limits<double>::infinity();
        double gamma_max = -std::numeric_limits<double>::infinity();
        double gamma_min = std::numeric_limits<double>::infinity();

        for (auto cell_index: cell.cells_) {
            const FVM::Entities::Cell_t& cell_k = mesh.cells_[cell_index];

            a_max = std::max(a_max, cell_k.a_hat_);
            a_min = std::min(a_min, cell_k.a_hat_);
            u_max = u_max.max(cell_k.u_hat_);
            u_min = u_min.min(cell_k.u_hat_);
            p_max = std::max(p_max, cell_k.p_hat_);
            p_min = std::min(p_min, cell_k.p_hat_);
            gamma_max = std::max(gamma_max, cell_k.gamma_hat_);
            gamma_min = std::min(gamma_min, cell_k.gamma_hat_);
        }

        for (auto face_index: cell.faces_) {
            const FVM::Entities::Face_t& face = mesh.faces_[face_index];

            const Vec2f delta = face.center_ - cell.center_;
            const double a_l = cell.a_hat_ + delta.dot(cell.a_derivative_);
            const Vec2f u_l = cell.u_hat_ + Vec2f(delta.dot(cell.ux_derivative_),
                                              delta.dot(cell.uy_derivative_));
            const double p_l = cell.p_hat_ + delta.dot(cell.p_derivative_);
            const double gamma_l = cell.gamma_hat_ + delta.dot(cell.gamma_derivative_);

            double y_a;
            double y_ux;
            double y_uy;
            double y_p;
            double y_gamma;

            if (a_l > cell.a_hat_) {
                y_a = (a_max - cell.a_hat_)/(a_l - cell.a_hat_);
            }
            else if (a_l < cell.a_hat_) {
                y_a = (a_min - cell.a_hat_)/(a_l - cell.a_hat_);
            }
            else {
                y_a = 2; // Will give phi = 1
            }

            if (u_l.x() > cell.u_hat_.x()) {
                y_ux = (u_max.x() - cell.u_hat_.x())/(u_l.x() - cell.u_hat_.x());
            }
            else if (u_l.x() < cell.u_hat_.x()) {
                y_ux = (u_min.x() - cell.u_hat_.x())/(u_l.x() - cell.u_hat_.x());
            }
            else {
                y_ux = 2; // Will give phi = 1
            }

            if (u_l.y() > cell.u_hat_.y()) {
                y_uy = (u_max.y() - cell.u_hat_.y())/(u_l.y() - cell.u_hat_.y());
            }
            else if (u_l.y() < cell.u_hat_.y()) {
                y_uy = (u_min.y() - cell.u_hat_.y())/(u_l.y() - cell.u_hat_.y());
            }
            else {
                y_uy = 2; // Will give phi = 1
            }

            if (p_l > cell.p_hat_) {
                y_p = (p_max - cell.p_hat_)/(p_l - cell.p_hat_);
            }
            else if (p_l < cell.p_hat_) {
                y_p = (p_min - cell.p_hat_)/(p_l - cell.p_hat_);
            }
            else {
                y_p = 2; // Will give phi = 1
            }

            if (gamma_l > cell.gamma_hat_) {
                y_gamma = (gamma_max - cell.gamma_hat_)/(gamma_l - cell.gamma_hat_);
            }
            else if (gamma_l < cell.gamma_hat_) {
                y_gamma = (gamma_min - cell.gamma_hat_)/(gamma_l - cell.gamma_hat_);
            }
            else {
                y_gamma = 2; // Will give phi = 1
            }

            const double phi_a_l = (std::pow(y_a, 2) + 2 * y_a)/(std::pow(y_a, 2) + y_a + 2);
            const Vec2f phi_u_l {(std::pow(y_ux, 2) + 2 * y_ux)/(std::pow(y_ux, 2) + y_ux + 2),
                                 (std::pow(y_uy, 2) + 2 * y_uy)/(std::pow(y_uy, 2) + y_uy + 2)};
            const double phi_p_l = (std::pow(y_p, 2) + 2 * y_p)/(std::pow(y_p, 2) + y_p + 2);
            const double phi_gamma_l = (std::pow(y_gamma, 2) + 2 * y_gamma)/(std::pow(y_gamma, 2) + y_gamma + 2);

            phi_a = std::min(phi_a, phi_a_l);
            phi_u.min(phi_u_l);
            phi_p = std::min(phi_p, phi_p_l);
            phi_gamma = std::min(phi_gamma, phi_gamma_l);
        }

        cell.a_derivative_ *= phi_a;
        cell.ux_derivative_ *= phi_u.x();
        cell.uy_derivative_ *= phi_u.y();
        cell.p_derivative_ *= phi_p;
        cell.gamma_derivative_ *= phi_gamma;
    }
}
