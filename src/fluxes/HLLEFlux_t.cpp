#include "fluxes/HLLEFlux_t.h"
#include <cmath>
#include <algorithm>

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
    
}
