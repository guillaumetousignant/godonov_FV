#include "fluxes/HLLEFlux_t.h"
#include <cmath>
#include <algorithm>

using FVM::Entities::Vec2f;

FVM::Fluxes::HLLEFlux_t::HLLEFlux_t() {}

FVM::Fluxes::HLLEFlux_t::~HLLEFlux_t() {}

void FVM::Fluxes::HLLEFlux_t::calculate_fluxes(double delta_t, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &F_1, std::vector<double> &F_2, std::vector<double> &F_3) {
    #pragma omp parallel for schedule(guided)
    for (int i = 0; i < F_1.size(); ++i) {
        // This would be all better if I stored those in the mesh, but that would not be drop-in replacable with the Riemann problem fluxes.
        const double F_1_L = gamma[i] * p[i] * u[i] /(std::pow(a[i], 2)); // Those should probably be cached somewhere, they are computed twice.
        const double F_2_L = gamma[i] * p[i] * std::pow(u[i], 2)/std::pow(a[i], 2) + p[i];
        const double F_3_L = u[i] * (gamma[i] * p[i] /(gamma[i] - 1) + gamma[i] * p[i] * std::pow(u[i], 2) * 0.5 /std::pow(a[i], 2));
        const double F_1_R = gamma[i+1] * p[i+1] * u[i+1] /(std::pow(a[i+1], 2));
        const double F_2_R = gamma[i+1] * p[i+1] * std::pow(u[i+1], 2)/std::pow(a[i+1], 2) + p[i+1];
        const double F_3_R = u[i+1] * (gamma[i+1] * p[i+1] /(gamma[i+1] - 1) + gamma[i+1] * p[i+1] * std::pow(u[i+1], 2) * 0.5 /std::pow(a[i+1], 2));

        const double U_1_L = gamma[i] * p[i] /std::pow(a[i], 2); // Those should probably be cached somewhere, they are computed twice.
        const double U_2_L = gamma[i] * p[i] * u[i]/std::pow(a[i], 2);
        const double U_3_L = p[i]/(gamma[i] - 1) + gamma[i] * p[i] * std::pow(u[i], 2) * 0.5/std::pow(a[i], 2);
        const double U_1_R = gamma[i+1] * p[i+1] /std::pow(a[i+1], 2);
        const double U_2_R = gamma[i+1] * p[i+1] * u[i+1]/std::pow(a[i+1], 2);
        const double U_3_R = p[i+1]/(gamma[i+1] - 1) + gamma[i+1] * p[i+1] * std::pow(u[i+1], 2) * 0.5/std::pow(a[i+1], 2);

        const double rho_p_L = std::sqrt(gamma[i] * p[i]); // We need it like 4 times. I assume it would have been optimized by the compiler anyway?
        const double rho_p_R = std::sqrt(gamma[i+1] * p[i+1]); // We need it like 4 times. I assume it would have been optimized by the compiler anyway?

        const double u_hat = ((rho_p_L * u[i]/a[i]) + (rho_p_R * u[i+1]/a[i+1])) / 
                                ((rho_p_L / a[i]) + (rho_p_R / a[i+1]));
        const double h_hat = (rho_p_L * (std::pow(a[i], 2) / (gamma[i] - 1) + std::pow(u[i], 2) * 0.5) /a[i] 
                                + rho_p_R * (std::pow(a[i+1], 2) / (gamma[i+1] - 1) + std::pow(u[i+1], 2) * 0.5) /a[i+1])
                                / ((rho_p_L / a[i]) + (rho_p_R / a[i+1]));
        //const double rho_hat = std::sqrt(gamma[i] * p[i] * gamma[i+1] * p[i+1])/(a[i] * a[i+1]);
        const double gamma_hat = (gamma[i] + gamma[i+1]) * 0.5; // Not sure, with an equation for this we could solve for a_hat directly.

        const double a_hat = std::sqrt((h_hat - std::pow(u_hat, 2) * 0.5) * (gamma_hat - 1));

        const double lambda_minus = std::min(u[i] - a[i], u_hat - a_hat); // Not sure about those
        const double lambda_plus = std::max(u[i+1] + a[i+1], u_hat + a_hat); // Not sure about those

        if (lambda_minus > 0) {
            F_1[i] = F_1_L;
            F_2[i] = F_2_L;
            F_3[i] = F_3_L;
        }
        else if (lambda_plus < 0) {
            F_1[i] = F_1_R;
            F_2[i] = F_2_R;
            F_3[i] = F_3_R;
        }
        else {
            F_1[i] = (lambda_plus * F_1_L - lambda_minus * F_1_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_1_R - U_1_L)/(lambda_plus - lambda_minus);
            F_2[i] = (lambda_plus * F_2_L - lambda_minus * F_2_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_2_R - U_2_L)/(lambda_plus - lambda_minus);
            F_3[i] = (lambda_plus * F_3_L - lambda_minus * F_3_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_3_R - U_3_L)/(lambda_plus - lambda_minus);
        }
    }
}

void FVM::Fluxes::HLLEFlux_t::calculate_fluxes_higher_order(double delta_t, const std::vector<double> x, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &F_1, std::vector<double> &F_2, std::vector<double> &F_3, const std::vector<double> du_dx, const std::vector<double> da_dx, const std::vector<double> dp_dx) {
    #pragma omp parallel for schedule(guided)
    for (int i = 0; i < F_1.size(); ++i) {
        const double delta_x = (x[i+1] - x[i])/2; // CHECK this won't work with non-uniform meshes! Use (x - x_i)

        const double gamma_L = gamma[i];
        const double u_L = u[i] + du_dx[i] * delta_x;
        const double a_L = a[i] + da_dx[i] * delta_x;
        const double p_L = p[i] + dp_dx[i] * delta_x;
        const double gamma_R = gamma[i+1];
        const double u_R = u[i+1] - du_dx[i+1] * delta_x;
        const double a_R = a[i+1] - da_dx[i+1] * delta_x;
        const double p_R = p[i+1] - dp_dx[i+1] * delta_x;

        // This would be all better if I stored those in the mesh, but that would not be drop-in replacable with the Riemann problem fluxes.
        const double F_1_L = gamma_L * p_L * u_L /(std::pow(a_L, 2)); // Those should probably be cached somewhere, they are computed twice.
        const double F_2_L = gamma_L * p_L * std::pow(u_L, 2)/std::pow(a_L, 2) + p_L;
        const double F_3_L = u_L * (gamma_L * p_L /(gamma_L - 1) + gamma_L * p_L * std::pow(u_L, 2) * 0.5 /std::pow(a_L, 2));
        const double F_1_R = gamma_R * p_R * u_R /(std::pow(a_R, 2));
        const double F_2_R = gamma_R * p_R * std::pow(u_R, 2)/std::pow(a_R, 2) + p_R;
        const double F_3_R = u_R * (gamma_R * p_R /(gamma_R - 1) + gamma_R * p_R * std::pow(u_R, 2) * 0.5 /std::pow(a_R, 2));

        const double U_1_L = gamma_L * p_L /std::pow(a_L, 2); // Those should probably be cached somewhere, they are computed twice.
        const double U_2_L = gamma_L * p_L * u_L/std::pow(a_L, 2);
        const double U_3_L = p_L/(gamma_L - 1) + gamma_L * p_L * std::pow(u_L, 2) * 0.5/std::pow(a_L, 2);
        const double U_1_R = gamma_R * p_R /std::pow(a_R, 2);
        const double U_2_R = gamma_R * p_R * u_R/std::pow(a_R, 2);
        const double U_3_R = p_R/(gamma_R - 1) + gamma_R * p_R * std::pow(u_R, 2) * 0.5/std::pow(a_R, 2);

        const double rho_p_L = std::sqrt(gamma_L * p_L); // We need it like 4 times. I assume it would have been optimized by the compiler anyway?
        const double rho_p_R = std::sqrt(gamma_R * p_R); // We need it like 4 times. I assume it would have been optimized by the compiler anyway?

        const double u_hat = ((rho_p_L * u_L/a_L) + (rho_p_R * u_R/a_R)) / 
                                ((rho_p_L / a_L) + (rho_p_R / a_R));
        const double h_hat = (rho_p_L * (std::pow(a_L, 2) / (gamma_L - 1) + std::pow(u_L, 2) * 0.5) /a_L 
                                + rho_p_R * (std::pow(a_R, 2) / (gamma_R - 1) + std::pow(u_R, 2) * 0.5) /a_R)
                                / ((rho_p_L / a_L) + (rho_p_R / a_R));
        //const double rho_hat = std::sqrt(gamma_L * p_L * gamma_R * p_R)/(a_L * a_R);
        const double gamma_hat = (gamma_L + gamma_R) * 0.5; // Not sure, with an equation for this we could solve for a_hat directly.

        const double a_hat = std::sqrt((h_hat - std::pow(u_hat, 2) * 0.5) * (gamma_hat - 1));

        const double lambda_minus = std::min(u_L - a_L, u_hat - a_hat); // Not sure about those
        const double lambda_plus = std::max(u_R + a_R, u_hat + a_hat); // Not sure about those

        if (lambda_minus > 0) {
            F_1[i] = F_1_L;
            F_2[i] = F_2_L;
            F_3[i] = F_3_L;
        }
        else if (lambda_plus < 0) {
            F_1[i] = F_1_R;
            F_2[i] = F_2_R;
            F_3[i] = F_3_R;
        }
        else {
            F_1[i] = (lambda_plus * F_1_L - lambda_minus * F_1_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_1_R - U_1_L)/(lambda_plus - lambda_minus);
            F_2[i] = (lambda_plus * F_2_L - lambda_minus * F_2_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_2_R - U_2_L)/(lambda_plus - lambda_minus);
            F_3[i] = (lambda_plus * F_3_L - lambda_minus * F_3_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_3_R - U_3_L)/(lambda_plus - lambda_minus);
        }
    }
}

void FVM::Fluxes::HLLEFlux_t::calculate_fluxes(double delta_t, FVM::Entities::Mesh2D_t &mesh) {
    #pragma omp parallel for schedule(guided)
    for (long long i = 0; i < mesh.faces_.size(); ++i) { // Because msvc openmp wants a signed index, so no size_t :(
        FVM::Entities::Face_t& face = mesh.faces_[i];
        const FVM::Entities::Cell_t& cell_L = mesh.cells_[face.cells_[0]];
        const FVM::Entities::Cell_t& cell_R = mesh.cells_[face.cells_[1]];

        const Vec2f u_prime_L(cell_L.u_.dot(face.normal_), cell_L.u_.dot(face.tangent_));
        const Vec2f u_prime_R(cell_R.u_.dot(face.normal_), cell_R.u_.dot(face.tangent_));

        const double gamma_p_L = std::sqrt(cell_L.gamma_ * cell_L.p_); // We need it like 4 times. I assume it would have been optimized by the compiler anyway?
        const double gamma_p_R = std::sqrt(cell_R.gamma_ * cell_R.p_); // We need it like 4 times. I assume it would have been optimized by the compiler anyway?

        const double u_hat = ((gamma_p_L * u_prime_L.x()/cell_L.a_) + (gamma_p_R * u_prime_R.x()/cell_R.a_)) / 
                                ((gamma_p_L / cell_L.a_) + (gamma_p_R / cell_R.a_));
        const double h_hat = (gamma_p_L * (std::pow(cell_L.a_, 2) / (cell_L.gamma_ - 1) + std::pow(u_prime_L.x(), 2) * 0.5) /cell_L.a_
                                + gamma_p_R * (std::pow(cell_R.a_, 2) / (cell_R.gamma_ - 1) + std::pow(u_prime_R.x(), 2) * 0.5) /cell_R.a_)
                                / ((gamma_p_L / cell_L.a_) + (gamma_p_R / cell_R.a_));
        //const double rho_hat = std::sqrt(gamma[i] * p[i] * gamma[i+1] * p[i+1])/(a[i] * a[i+1]);
        const double gamma_hat = (cell_L.gamma_ + cell_R.gamma_) * 0.5; // Not sure, with an equation for this we could solve for a_hat directly.

        const double a_hat = std::sqrt((h_hat - std::pow(u_hat, 2) * 0.5) * (gamma_hat - 1));

        const double lambda_minus = std::min(u_prime_L.x() - cell_L.a_, u_hat - a_hat); // Not sure about those
        const double lambda_plus = std::max(u_prime_R.x() + cell_R.a_, u_hat + a_hat); // Not sure about those

        const Vec2f normal_inv = Vec2f(face.normal_[0], face.tangent_[0]);
        const Vec2f tangent_inv = Vec2f(face.normal_[1], face.tangent_[1]);

        const double U_1_L = cell_L.gamma_ * cell_L.p_ /std::pow(cell_L.a_, 2); // Those should probably be cached somewhere, they are computed twice.
        const double U_2_L = cell_L.gamma_ * cell_L.p_ * u_prime_L.x()/std::pow(cell_L.a_, 2);
        const double U_3_L = cell_L.gamma_ * cell_L.p_ * u_prime_L.y()/std::pow(cell_L.a_, 2);
        const double U_4_L = cell_L.p_/(cell_L.gamma_ - 1) + cell_L.gamma_ * cell_L.p_ * u_prime_L.magnitudeSquared() * 0.5/std::pow(cell_L.a_, 2);
        const double U_1_R = cell_R.gamma_ * cell_R.p_ /std::pow(cell_R.a_, 2); // Those should probably be cached somewhere, they are computed twice.
        const double U_2_R = cell_R.gamma_ * cell_R.p_ * u_prime_R.x()/std::pow(cell_R.a_, 2);
        const double U_3_R = cell_R.gamma_ * cell_R.p_ * u_prime_R.y()/std::pow(cell_R.a_, 2);
        const double U_4_R = cell_R.p_/(cell_R.gamma_ - 1) + cell_R.gamma_ * cell_R.p_ * u_prime_R.magnitudeSquared() * 0.5/std::pow(cell_R.a_, 2);

        // This would be all better if I stored those in the mesh, but that would not be drop-in replacable with the Riemann problem fluxes.
        const double F_1_L = u_prime_L.x() * cell_L.gamma_ * cell_L.p_ /(std::pow(cell_L.a_, 2)); // Those should probably be cached somewhere, they are computed twice.
        const double F_2_L = cell_L.gamma_ * cell_L.p_ * std::pow(u_prime_L.x(), 2)/std::pow(cell_L.a_, 2) + cell_L.p_;
        const double F_3_L = u_prime_L.x() * u_prime_L.y() * cell_L.gamma_ * cell_L.p_/std::pow(cell_L.a_, 2);
        const double F_4_L = u_prime_L.x() * (cell_L.gamma_ * cell_L.p_ /(cell_L.gamma_ - 1) + cell_L.gamma_ * cell_L.p_ * (std::pow(u_prime_L.x(), 2) + std::pow(u_prime_L.y(), 2)) * 0.5 /std::pow(cell_L.a_, 2));
        const double F_1_R = u_prime_R.x() * cell_R.gamma_ * cell_R.p_ /(std::pow(cell_R.a_, 2)); // Those should probably be cached somewhere, they are computed twice.
        const double F_2_R = cell_R.gamma_ * cell_R.p_ * std::pow(u_prime_R.x(), 2)/std::pow(cell_R.a_, 2) + cell_R.p_;
        const double F_3_R = u_prime_R.x() * u_prime_R.y() * cell_R.gamma_ * cell_R.p_/std::pow(cell_R.a_, 2);
        const double F_4_R = u_prime_R.x() * (cell_R.gamma_ * cell_R.p_ /(cell_R.gamma_ - 1) + cell_R.gamma_ * cell_R.p_ * (std::pow(u_prime_R.x(), 2) + std::pow(u_prime_R.y(), 2)) * 0.5 /std::pow(cell_R.a_, 2));


        if (lambda_minus > 0) {
            face.F_1_ = F_1_L;
            face.F_2_ = normal_inv.dot({F_2_L, F_3_L}); //
            face.F_3_ = tangent_inv.dot({F_2_L, F_3_L});
            face.F_4_ = F_4_L;
        }
        else if (lambda_plus < 0) {
            face.F_1_ = F_1_R;
            face.F_2_ = normal_inv.dot({F_2_R, F_3_R});
            face.F_3_ = tangent_inv.dot({F_2_R, F_3_R});
            face.F_4_ = F_4_R;
        }
        else {
            const double F_2 = (lambda_plus * F_2_L - lambda_minus * F_2_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_2_R - U_2_L)/(lambda_plus - lambda_minus);
            const double F_3 = (lambda_plus * F_3_L - lambda_minus * F_3_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_3_R - U_3_L)/(lambda_plus - lambda_minus);

            face.F_1_ = (lambda_plus * F_1_L - lambda_minus * F_1_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_1_R - U_1_L)/(lambda_plus - lambda_minus);
            face.F_2_ = normal_inv.dot({F_2, F_3});
            face.F_3_ = tangent_inv.dot({F_2, F_3});
            face.F_4_ = (lambda_plus * F_4_L - lambda_minus * F_4_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_4_R - U_4_L)/(lambda_plus - lambda_minus);
        }
    }
}

void FVM::Fluxes::HLLEFlux_t::calculate_fluxes_higher_order(double delta_t, FVM::Entities::Mesh2D_t &mesh) {
    #pragma omp parallel for schedule(guided)
    for (long long i = 0; i < mesh.faces_.size(); ++i) { // Because msvc openmp wants a signed index, so no size_t :(
        FVM::Entities::Face_t& face = mesh.faces_[i];
        const FVM::Entities::Cell_t& cell_L = mesh.cells_[face.cells_[0]];
        const FVM::Entities::Cell_t& cell_R = mesh.cells_[face.cells_[1]];

        const Vec2f delta_L = face.center_ - cell_L.center_;
        const Vec2f delta_R = face.center_ - cell_R.center_;

        const double gamma_L = cell_L.gamma_ + delta_L.dot(cell_L.gamma_derivative_);
        const Vec2f u_L = cell_L.u_ + Vec2f(delta_L.dot(cell_L.ux_derivative_), delta_L.dot(cell_L.uy_derivative_));
        const double a_L = cell_L.a_ + delta_L.dot(cell_L.a_derivative_);
        const double p_L = cell_L.p_ + delta_L.dot(cell_L.p_derivative_);
        const double gamma_R = cell_R.gamma_ + delta_R.dot(cell_R.gamma_derivative_);
        const Vec2f u_R = cell_R.u_ + Vec2f(delta_R.dot(cell_R.ux_derivative_), delta_R.dot(cell_R.uy_derivative_));
        const double a_R = cell_R.a_ + delta_R.dot(cell_R.a_derivative_);
        const double p_R = cell_R.p_ + delta_R.dot(cell_R.p_derivative_);

        const Vec2f u_prime_L(u_L.dot(face.normal_), u_L.dot(face.tangent_));
        const Vec2f u_prime_R(u_R.dot(face.normal_), u_R.dot(face.tangent_));

        const double gamma_p_L = std::sqrt(gamma_L * p_L); // We need it like 4 times. I assume it would have been optimized by the compiler anyway?
        const double gamma_p_R = std::sqrt(gamma_R * p_R); // We need it like 4 times. I assume it would have been optimized by the compiler anyway?

        const double u_hat = ((gamma_p_L * u_prime_L.x()/a_L) + (gamma_p_R * u_prime_R.x()/a_R)) / 
                                ((gamma_p_L / a_L) + (gamma_p_R / a_R));
        const double h_hat = (gamma_p_L * (std::pow(a_L, 2) / (gamma_L - 1) + std::pow(u_prime_L.x(), 2) * 0.5) /a_L
                                + gamma_p_R * (std::pow(a_R, 2) / (gamma_R - 1) + std::pow(u_prime_R.x(), 2) * 0.5) /a_R)
                                / ((gamma_p_L / a_L) + (gamma_p_R / a_R));
        //const double rho_hat = std::sqrt(gamma[i] * p[i] * gamma[i+1] * p[i+1])/(a[i] * a[i+1]);
        const double gamma_hat = (gamma_L + gamma_R) * 0.5; // Not sure, with an equation for this we could solve for a_hat directly.

        const double a_hat = std::sqrt((h_hat - std::pow(u_hat, 2) * 0.5) * (gamma_hat - 1));

        const double lambda_minus = std::min(u_prime_L.x() - a_L, u_hat - a_hat); // Not sure about those
        const double lambda_plus = std::max(u_prime_R.x() + a_R, u_hat + a_hat); // Not sure about those

        const Vec2f normal_inv = Vec2f(face.normal_[0], face.tangent_[0]);
        const Vec2f tangent_inv = Vec2f(face.normal_[1], face.tangent_[1]);

        const double U_1_L = gamma_L * p_L /std::pow(a_L, 2); // Those should probably be cached somewhere, they are computed twice.
        const double U_2_L = gamma_L * p_L * u_prime_L.x()/std::pow(a_L, 2);
        const double U_3_L = gamma_L * p_L * u_prime_L.y()/std::pow(a_L, 2);
        const double U_4_L = p_L/(gamma_L - 1) + gamma_L * p_L * u_prime_L.magnitudeSquared() * 0.5/std::pow(a_L, 2);
        const double U_1_R = gamma_R * p_R /std::pow(a_R, 2); // Those should probably be cached somewhere, they are computed twice.
        const double U_2_R = gamma_R * p_R * u_prime_R.x()/std::pow(a_R, 2);
        const double U_3_R = gamma_R * p_R * u_prime_R.y()/std::pow(a_R, 2);
        const double U_4_R = p_R/(gamma_R - 1) + gamma_R * p_R * u_prime_R.magnitudeSquared() * 0.5/std::pow(a_R, 2);

        // This would be all better if I stored those in the mesh, but that would not be drop-in replacable with the Riemann problem fluxes.
        const double F_1_L = u_prime_L.x() * gamma_L * p_L /(std::pow(a_L, 2)); // Those should probably be cached somewhere, they are computed twice.
        const double F_2_L = gamma_L * p_L * std::pow(u_prime_L.x(), 2)/std::pow(a_L, 2) + p_L;
        const double F_3_L = u_prime_L.x() * u_prime_L.y() * gamma_L * p_L/std::pow(a_L, 2);
        const double F_4_L = u_prime_L.x() * (gamma_L * p_L /(gamma_L - 1) + gamma_L * p_L * (std::pow(u_prime_L.x(), 2) + std::pow(u_prime_L.y(), 2)) * 0.5 /std::pow(a_L, 2));
        const double F_1_R = u_prime_R.x() * gamma_R * p_R /(std::pow(a_R, 2)); // Those should probably be cached somewhere, they are computed twice.
        const double F_2_R = gamma_R * p_R * std::pow(u_prime_R.x(), 2)/std::pow(a_R, 2) + p_R;
        const double F_3_R = u_prime_R.x() * u_prime_R.y() * gamma_R * p_R/std::pow(a_R, 2);
        const double F_4_R = u_prime_R.x() * (gamma_R * p_R /(gamma_R - 1) + gamma_R * p_R * (std::pow(u_prime_R.x(), 2) + std::pow(u_prime_R.y(), 2)) * 0.5 /std::pow(a_R, 2));


        if (lambda_minus > 0) {
            face.F_1_ = F_1_L;
            face.F_2_ = normal_inv.dot({F_2_L, F_3_L}); //
            face.F_3_ = tangent_inv.dot({F_2_L, F_3_L});
            face.F_4_ = F_4_L;
        }
        else if (lambda_plus < 0) {
            face.F_1_ = F_1_R;
            face.F_2_ = normal_inv.dot({F_2_R, F_3_R});
            face.F_3_ = tangent_inv.dot({F_2_R, F_3_R});
            face.F_4_ = F_4_R;
        }
        else {
            const double F_2 = (lambda_plus * F_2_L - lambda_minus * F_2_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_2_R - U_2_L)/(lambda_plus - lambda_minus);
            const double F_3 = (lambda_plus * F_3_L - lambda_minus * F_3_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_3_R - U_3_L)/(lambda_plus - lambda_minus);

            face.F_1_ = (lambda_plus * F_1_L - lambda_minus * F_1_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_1_R - U_1_L)/(lambda_plus - lambda_minus);
            face.F_2_ = normal_inv.dot({F_2, F_3});
            face.F_3_ = tangent_inv.dot({F_2, F_3});
            face.F_4_ = (lambda_plus * F_4_L - lambda_minus * F_4_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_4_R - U_4_L)/(lambda_plus - lambda_minus);
        }
    }
}

void FVM::Fluxes::HLLEFlux_t::calculate_fluxes_higher_order_hat(double delta_t, FVM::Entities::Mesh2D_t &mesh) {
    #pragma omp parallel for schedule(guided)
    for (long long i = 0; i < mesh.faces_.size(); ++i) { // Because msvc openmp wants a signed index, so no size_t :(
        FVM::Entities::Face_t& face = mesh.faces_[i];
        const FVM::Entities::Cell_t& cell_L = mesh.cells_[face.cells_[0]];
        const FVM::Entities::Cell_t& cell_R = mesh.cells_[face.cells_[1]];

        const Vec2f delta_L = face.center_ - cell_L.center_;
        const Vec2f delta_R = face.center_ - cell_R.center_;

        const double gamma_L = cell_L.gamma_hat_ + delta_L.dot(cell_L.gamma_derivative_);
        const Vec2f u_L = cell_L.u_hat_ + Vec2f(delta_L.dot(cell_L.ux_derivative_), delta_L.dot(cell_L.uy_derivative_));
        const double a_L = cell_L.a_hat_ + delta_L.dot(cell_L.a_derivative_);
        const double p_L = cell_L.p_hat_ + delta_L.dot(cell_L.p_derivative_);
        const double gamma_R = cell_R.gamma_hat_ + delta_R.dot(cell_R.gamma_derivative_);
        const Vec2f u_R = cell_R.u_hat_ + Vec2f(delta_R.dot(cell_R.ux_derivative_), delta_R.dot(cell_R.uy_derivative_));
        const double a_R = cell_R.a_hat_ + delta_R.dot(cell_R.a_derivative_);
        const double p_R = cell_R.p_hat_ + delta_R.dot(cell_R.p_derivative_);

        const Vec2f u_prime_L(u_L.dot(face.normal_), u_L.dot(face.tangent_));
        const Vec2f u_prime_R(u_R.dot(face.normal_), u_R.dot(face.tangent_));

        const double gamma_p_L = std::sqrt(gamma_L * p_L); // We need it like 4 times. I assume it would have been optimized by the compiler anyway?
        const double gamma_p_R = std::sqrt(gamma_R * p_R); // We need it like 4 times. I assume it would have been optimized by the compiler anyway?

        const double u_hat = ((gamma_p_L * u_prime_L.x()/a_L) + (gamma_p_R * u_prime_R.x()/a_R)) / 
                                ((gamma_p_L / a_L) + (gamma_p_R / a_R));
        const double h_hat = (gamma_p_L * (std::pow(a_L, 2) / (gamma_L - 1) + std::pow(u_prime_L.x(), 2) * 0.5) /a_L
                                + gamma_p_R * (std::pow(a_R, 2) / (gamma_R - 1) + std::pow(u_prime_R.x(), 2) * 0.5) /a_R)
                                / ((gamma_p_L / a_L) + (gamma_p_R / a_R));
        //const double rho_hat = std::sqrt(gamma[i] * p[i] * gamma[i+1] * p[i+1])/(a[i] * a[i+1]);
        const double gamma_hat = (gamma_L + gamma_R) * 0.5; // Not sure, with an equation for this we could solve for a_hat directly.

        const double a_hat = std::sqrt((h_hat - std::pow(u_hat, 2) * 0.5) * (gamma_hat - 1));

        const double lambda_minus = std::min(u_prime_L.x() - a_L, u_hat - a_hat); // Not sure about those
        const double lambda_plus = std::max(u_prime_R.x() + a_R, u_hat + a_hat); // Not sure about those

        const Vec2f normal_inv = Vec2f(face.normal_[0], face.tangent_[0]);
        const Vec2f tangent_inv = Vec2f(face.normal_[1], face.tangent_[1]);

        const double U_1_L = gamma_L * p_L /std::pow(a_L, 2); // Those should probably be cached somewhere, they are computed twice.
        const double U_2_L = gamma_L * p_L * u_prime_L.x()/std::pow(a_L, 2);
        const double U_3_L = gamma_L * p_L * u_prime_L.y()/std::pow(a_L, 2);
        const double U_4_L = p_L/(gamma_L - 1) + gamma_L * p_L * u_prime_L.magnitudeSquared() * 0.5/std::pow(a_L, 2);
        const double U_1_R = gamma_R * p_R /std::pow(a_R, 2); // Those should probably be cached somewhere, they are computed twice.
        const double U_2_R = gamma_R * p_R * u_prime_R.x()/std::pow(a_R, 2);
        const double U_3_R = gamma_R * p_R * u_prime_R.y()/std::pow(a_R, 2);
        const double U_4_R = p_R/(gamma_R - 1) + gamma_R * p_R * u_prime_R.magnitudeSquared() * 0.5/std::pow(a_R, 2);

        // This would be all better if I stored those in the mesh, but that would not be drop-in replacable with the Riemann problem fluxes.
        const double F_1_L = u_prime_L.x() * gamma_L * p_L /(std::pow(a_L, 2)); // Those should probably be cached somewhere, they are computed twice.
        const double F_2_L = gamma_L * p_L * std::pow(u_prime_L.x(), 2)/std::pow(a_L, 2) + p_L;
        const double F_3_L = u_prime_L.x() * u_prime_L.y() * gamma_L * p_L/std::pow(a_L, 2);
        const double F_4_L = u_prime_L.x() * (gamma_L * p_L /(gamma_L - 1) + gamma_L * p_L * (std::pow(u_prime_L.x(), 2) + std::pow(u_prime_L.y(), 2)) * 0.5 /std::pow(a_L, 2));
        const double F_1_R = u_prime_R.x() * gamma_R * p_R /(std::pow(a_R, 2)); // Those should probably be cached somewhere, they are computed twice.
        const double F_2_R = gamma_R * p_R * std::pow(u_prime_R.x(), 2)/std::pow(a_R, 2) + p_R;
        const double F_3_R = u_prime_R.x() * u_prime_R.y() * gamma_R * p_R/std::pow(a_R, 2);
        const double F_4_R = u_prime_R.x() * (gamma_R * p_R /(gamma_R - 1) + gamma_R * p_R * (std::pow(u_prime_R.x(), 2) + std::pow(u_prime_R.y(), 2)) * 0.5 /std::pow(a_R, 2));


        if (lambda_minus > 0) {
            face.F_1_hat_ = F_1_L;
            face.F_2_hat_ = normal_inv.dot({F_2_L, F_3_L}); //
            face.F_3_hat_ = tangent_inv.dot({F_2_L, F_3_L});
            face.F_4_hat_ = F_4_L;
        }
        else if (lambda_plus < 0) {
            face.F_1_hat_ = F_1_R;
            face.F_2_hat_ = normal_inv.dot({F_2_R, F_3_R});
            face.F_3_hat_ = tangent_inv.dot({F_2_R, F_3_R});
            face.F_4_hat_ = F_4_R;
        }
        else {
            const double F_2 = (lambda_plus * F_2_L - lambda_minus * F_2_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_2_R - U_2_L)/(lambda_plus - lambda_minus);
            const double F_3 = (lambda_plus * F_3_L - lambda_minus * F_3_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_3_R - U_3_L)/(lambda_plus - lambda_minus);

            face.F_1_hat_ = (lambda_plus * F_1_L - lambda_minus * F_1_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_1_R - U_1_L)/(lambda_plus - lambda_minus);
            face.F_2_hat_ = normal_inv.dot({F_2, F_3});
            face.F_3_hat_ = tangent_inv.dot({F_2, F_3});
            face.F_4_hat_ = (lambda_plus * F_4_L - lambda_minus * F_4_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_4_R - U_4_L)/(lambda_plus - lambda_minus);
        }
    }
}
