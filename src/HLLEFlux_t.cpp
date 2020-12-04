#include "HLLEFlux_t.h"
#include <cmath>
#include <algorithm>

HLLEFlux_t::HLLEFlux_t(int n_faces) {}

HLLEFlux_t::~HLLEFlux_t() {}

void HLLEFlux_t::calculate_fluxes(Mesh1D_t &mesh, double delta_t) {
    for (int i = 0; i < mesh.n_faces_; ++i) {
        // This would be all better if I stored those in the mesh, but that would not be drop-in replacable with the Riemann problem fluxes.
        const double F_1_L = mesh.gamma_[i] * mesh.p_[i] * mesh.u_[i] /(std::pow(mesh.a_[i], 2)); // Those should probably be cached somewhere, they are computed twice.
        const double F_2_L = mesh.gamma_[i] * mesh.p_[i] * std::pow(mesh.u_[i], 2)/std::pow(mesh.a_[i], 2) + mesh.p_[i];
        const double F_3_L = mesh.u_[i] * (mesh.gamma_[i] * mesh.p_[i] /(mesh.gamma_[i] - 1) + mesh.gamma_[i] * mesh.p_[i] * std::pow(mesh.u_[i], 2) * 0.5 /std::pow(mesh.a_[i], 2));
        const double F_1_R = mesh.gamma_[i+1] * mesh.p_[i+1] * mesh.u_[i+1] /(std::pow(mesh.a_[i+1], 2));
        const double F_2_R = mesh.gamma_[i+1] * mesh.p_[i+1] * std::pow(mesh.u_[i+1], 2)/std::pow(mesh.a_[i+1], 2) + mesh.p_[i+1];
        const double F_3_R = mesh.u_[i+1] * (mesh.gamma_[i+1] * mesh.p_[i+1] /(mesh.gamma_[i+1] - 1) + mesh.gamma_[i+1] * mesh.p_[i+1] * std::pow(mesh.u_[i+1], 2) * 0.5 /std::pow(mesh.a_[i+1], 2));

        const double U_1_L = mesh.gamma_[i] * mesh.p_[i] /std::pow(mesh.a_[i], 2); // Those should probably be cached somewhere, they are computed twice.
        const double U_2_L = mesh.gamma_[i] * mesh.p_[i] * mesh.u_[i]/std::pow(mesh.a_[i], 2);
        const double U_3_L = mesh.p_[i]/(mesh.gamma_[i] - 1) + mesh.gamma_[i] * mesh.p_[i] * std::pow(mesh.u_[i], 2) * 0.5/std::pow(mesh.a_[i], 2);
        const double U_1_R = mesh.gamma_[i+1] * mesh.p_[i+1] /std::pow(mesh.a_[i+1], 2);
        const double U_2_R = mesh.gamma_[i+1] * mesh.p_[i+1] * mesh.u_[i+1]/std::pow(mesh.a_[i+1], 2);
        const double U_3_R = mesh.p_[i+1]/(mesh.gamma_[i+1] - 1) + mesh.gamma_[i+1] * mesh.p_[i+1] * std::pow(mesh.u_[i+1], 2) * 0.5/std::pow(mesh.a_[i+1], 2);

        const double rho_p_L = std::sqrt(mesh.gamma_[i] * mesh.p_[i]); // We need it like 4 times. I assume it would have been optimized by the compiler anyway?
        const double rho_p_R = std::sqrt(mesh.gamma_[i+1] * mesh.p_[i+1]); // We need it like 4 times. I assume it would have been optimized by the compiler anyway?

        const double u_hat = ((rho_p_L * mesh.u_[i]/mesh.a_[i]) + (rho_p_R * mesh.u_[i+1]/mesh.a_[i+1])) / 
                                ((rho_p_L / mesh.a_[i]) + (rho_p_R / mesh.a_[i+1]));
        const double h_hat = (rho_p_L * (std::pow(mesh.a_[i], 2) / (mesh.gamma_[i] - 1) + std::pow(mesh.u_[i], 2) * 0.5) /mesh.a_[i] 
                                + rho_p_R * (std::pow(mesh.a_[i+1], 2) / (mesh.gamma_[i+1] - 1) + std::pow(mesh.u_[i+1], 2) * 0.5) /mesh.a_[i+1])
                                / ((rho_p_L / mesh.a_[i]) + (rho_p_R / mesh.a_[i+1]));
        //const double rho_hat = std::sqrt(mesh.gamma_[i] * mesh.p_[i] * mesh.gamma_[i+1] * mesh.p_[i+1])/(mesh.a_[i] * mesh.a_[i+1]);
        const double gamma_hat = (mesh.gamma_[i] + mesh.gamma_[i+1]) * 0.5; // Not sure, with an equation for this we could solve for a_hat directly.

        const double a_hat = std::sqrt((h_hat - std::pow(u_hat, 2) * 0.5) * (gamma_hat - 1));

        const double lambda_minus = std::min(mesh.u_[i] - mesh.a_[i], u_hat - a_hat); // Not sure about those
        const double lambda_plus = std::max(mesh.u_[i+1] + mesh.a_[i+1], u_hat + a_hat); // Not sure about those

        if (lambda_minus > 0) {
            mesh.F_1_[i] = F_1_L;
            mesh.F_2_[i] = F_2_L;
            mesh.F_3_[i] = F_3_L;
        }
        else if (lambda_plus < 0) {
            mesh.F_1_[i] = F_1_R;
            mesh.F_2_[i] = F_2_R;
            mesh.F_3_[i] = F_3_R;
        }
        else {
            mesh.F_1_[i] = (lambda_plus * F_1_L - lambda_minus * F_1_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_1_R - U_1_L)/(lambda_plus - lambda_minus);
            mesh.F_2_[i] = (lambda_plus * F_2_L - lambda_minus * F_2_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_2_R - U_2_L)/(lambda_plus - lambda_minus);
            mesh.F_3_[i] = (lambda_plus * F_3_L - lambda_minus * F_3_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_3_R - U_3_L)/(lambda_plus - lambda_minus);
        }
    }
}