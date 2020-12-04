#include "RoeEntropyFlux_t.h"

RoeEntropyFlux_t::RoeEntropyFlux_t(int n_faces) {}

RoeEntropyFlux_t::~RoeEntropyFlux_t() {}

void RoeEntropyFlux_t::invert_matrix(const double (&input)[9], double (&output)[9]) {
    const double invdet = 1/(input[0] * (input[4] * input[8] - input[7] * input[5]) -
                             input[1] * (input[3] * input[8] - input[5] * input[6]) +
                             input[2] * (input[3] * input[7] - input[4] * input[6]));

    output[0] = (input[4] * input[8] - input[7] * input[5]) * invdet;
    output[1] = (input[2] * input[7] - input[1] * input[8]) * invdet;
    output[2] = (input[1] * input[5] - input[2] * input[4]) * invdet;
    output[3] = (input[5] * input[6] - input[3] * input[8]) * invdet;
    output[4] = (input[0] * input[8] - input[2] * input[6]) * invdet;
    output[5] = (input[3] * input[2] - input[0] * input[3]) * invdet;
    output[6] = (input[3] * input[7] - input[6] * input[4]) * invdet;
    output[7] = (input[6] * input[1] - input[0] * input[7]) * invdet;
    output[8] = (input[0] * input[4] - input[3] * input[1]) * invdet;
}

void RoeEntropyFlux_t::multiply_matrix(const double (&left)[9], const double (&right)[9], double (&result)[9]) { // What do you mean loops
    result[0] = left[0] * right[0] + left[1] * right[3] + left[2] * right[6];
    result[1] = left[0] * right[1] + left[1] * right[4] + left[2] * right[7];
    result[2] = left[0] * right[2] + left[1] * right[5] + left[2] * right[8];
    result[3] = left[3] * right[0] + left[4] * right[3] + left[5] * right[6];
    result[4] = left[3] * right[1] + left[4] * right[4] + left[5] * right[7];
    result[5] = left[3] * right[2] + left[4] * right[5] + left[5] * right[8];
    result[6] = left[6] * right[0] + left[7] * right[3] + left[8] * right[6];
    result[7] = left[6] * right[1] + left[7] * right[4] + left[8] * right[7];
    result[8] = left[6] * right[2] + left[7] * right[5] + left[8] * right[8];
}


void RoeEntropyFlux_t::calculate_fluxes(Mesh1D_t &mesh, double delta_t) {
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

        double R_hat[9] = {1.0,                   1.0,                      1.0,
                           u_hat - a_hat,         u_hat,                    u_hat + a_hat,
                           h_hat - u_hat * a_hat, std::pow(u_hat, 2) * 0.5, h_hat + u_hat * a_hat};
        double Lambda_hat[9] = {std::abs(u_hat - a_hat), 0.0,             0.0,
                                      0.0,               std::abs(u_hat), 0.0,
                                      0.0,               0.0,             std::abs(u_hat + a_hat)};
        double L_hat[9];
        invert_matrix(R_hat, L_hat);

        double mult_matrix[9]; // Temp value, Roe matrix will be stored in R_hat
        multiply_matrix(R_hat, Lambda_hat, mult_matrix);
        multiply_matrix(mult_matrix, L_hat, R_hat); // Now R_hat is Roe matrix

        mesh.F_1_[i] = 0.5 * (F_1_R + F_1_L) - 0.5 * (R_hat[0] * (U_1_R - U_1_L) + R_hat[1] * (U_2_R - U_2_L) + R_hat[2] * (U_3_R - U_3_L));
        mesh.F_2_[i] = 0.5 * (F_2_R + F_2_L) - 0.5 * (R_hat[3] * (U_1_R - U_1_L) + R_hat[4] * (U_2_R - U_2_L) + R_hat[5] * (U_3_R - U_3_L));
        mesh.F_3_[i] = 0.5 * (F_3_R + F_3_L) - 0.5 * (R_hat[6] * (U_1_R - U_1_L) + R_hat[7] * (U_2_R - U_2_L) + R_hat[8] * (U_3_R - U_3_L));
    }
}