#include "RoeFlux_t.h"
#include <cmath>

RoeFlux_t::RoeFlux_t(int n_faces) {}

RoeFlux_t::~RoeFlux_t() {}

void RoeFlux_t::calculate_fluxes(Mesh1D_t &mesh, double delta_t) {
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


        const double lambda_minus = mesh.u_[i] - mesh.a_[i]; // Not sure about those
        const double lambda_plus = mesh.u_[i+1] + mesh.a_[i+1]; // Not sure about those

        mesh.F_1_[i] = (lambda_plus * F_1_L - lambda_minus * F_1_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_1_R - U_1_L)/(lambda_plus - lambda_minus);
        mesh.F_2_[i] = (lambda_plus * F_2_L - lambda_minus * F_2_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_2_R - U_2_L)/(lambda_plus - lambda_minus);
        mesh.F_3_[i] = (lambda_plus * F_3_L - lambda_minus * F_3_R)/(lambda_plus - lambda_minus) + lambda_plus * lambda_minus * (U_3_R - U_3_L)/(lambda_plus - lambda_minus);
    }
}