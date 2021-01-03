#include "fluxes/RoeFlux_t.h"
#include <cmath>

FVM::Fluxes::RoeFlux_t::RoeFlux_t() {}

FVM::Fluxes::RoeFlux_t::~RoeFlux_t() {}

void FVM::Fluxes::RoeFlux_t::invert_matrix(const double (&input)[9], double (&output)[9]) {
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

void FVM::Fluxes::RoeFlux_t::multiply_matrix(const double (&left)[9], const double (&right)[9], double (&result)[9]) { // What do you mean loops
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

void FVM::Fluxes::RoeFlux_t::calculate_fluxes(double delta_t, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &F_1, std::vector<double> &F_2, std::vector<double> &F_3) {
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

        F_1[i] = 0.5 * (F_1_R + F_1_L) - 0.5 * (R_hat[0] * (U_1_R - U_1_L) + R_hat[1] * (U_2_R - U_2_L) + R_hat[2] * (U_3_R - U_3_L));
        F_2[i] = 0.5 * (F_2_R + F_2_L) - 0.5 * (R_hat[3] * (U_1_R - U_1_L) + R_hat[4] * (U_2_R - U_2_L) + R_hat[5] * (U_3_R - U_3_L));
        F_3[i] = 0.5 * (F_3_R + F_3_L) - 0.5 * (R_hat[6] * (U_1_R - U_1_L) + R_hat[7] * (U_2_R - U_2_L) + R_hat[8] * (U_3_R - U_3_L));
    }
}

void FVM::Fluxes::RoeFlux_t::calculate_fluxes_higher_order(double delta_t, const std::vector<double> x, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &F_1, std::vector<double> &F_2, std::vector<double> &F_3, const std::vector<double> du_dx, const std::vector<double> da_dx, const std::vector<double> dp_dx) {
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

        F_1[i] = 0.5 * (F_1_R + F_1_L) - 0.5 * (R_hat[0] * (U_1_R - U_1_L) + R_hat[1] * (U_2_R - U_2_L) + R_hat[2] * (U_3_R - U_3_L));
        F_2[i] = 0.5 * (F_2_R + F_2_L) - 0.5 * (R_hat[3] * (U_1_R - U_1_L) + R_hat[4] * (U_2_R - U_2_L) + R_hat[5] * (U_3_R - U_3_L));
        F_3[i] = 0.5 * (F_3_R + F_3_L) - 0.5 * (R_hat[6] * (U_1_R - U_1_L) + R_hat[7] * (U_2_R - U_2_L) + R_hat[8] * (U_3_R - U_3_L));
    }
}

void FVM::Fluxes::RoeFlux_t::calculate_fluxes(double delta_t, FVM::Entities::Mesh2D_t &mesh) {
    
}

void FVM::Fluxes::RoeFlux_t::calculate_fluxes_higher_order(double delta_t, FVM::Entities::Mesh2D_t &mesh) {

}

void FVM::Fluxes::RoeFlux_t::calculate_fluxes_higher_order_hat(double delta_t, FVM::Entities::Mesh2D_t &mesh) {

}
