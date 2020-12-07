#ifndef FVM_RIEMANNPROBLEM_H
#define FVM_RIEMANNPROBLEM_H

namespace RiemannProblem {
    double u_star_initial_guess(double a_L, double a_R, double u_L, double u_R, double p_L, double p_R, double gamma_L, double gamma_R);
    double calculate_C(double a, double p, double gamma);
    void calculate_a_star(double a_L, double a_R, double u_L, double u_R, double p_L, double p_R, double gamma_L, double gamma_R, double u_star, double &a_star_L, double &a_star_R, double p_star_L, double p_star_R);
    void left_shock(double a_L, double u_L, double p_L, double gamma_L, double C_L, double u_star, double &p_star_L, double &p_star_prime_L);
    void right_shock(double a_R, double u_R, double p_R, double gamma_R, double C_R, double u_star, double &p_star_R, double &p_star_prime_R);
    void left_rarefaction(double a_L, double u_L, double p_L, double gamma_L, double u_star, double &a_star_L, double &p_star_L, double &p_star_prime_L);
    void right_rarefaction(double a_R, double u_R, double p_R, double gamma_R, double u_star, double &a_star_R, double &p_star_R, double &p_star_prime_R);
    void solve_flux(double a_L, double a_R, double u_L, double u_R, double p_L, double p_R, double gamma_L, double gamma_R, double C_L, double C_R, double &u_star, double &a_star_L, double &a_star_R, double &p_star_L, double &p_star_R, double &p_star_prime_L, double &p_star_prime_R);
    void get_boundary_state(double &a_face, double &u_face, double &p_face, double &gamma_face, double delta_t, double a_L, double a_R, double u_L, double u_R, double p_L, double p_R, double gamma_L, double gamma_R, double u_star, double a_star_L, double a_star_R, double p_star_L, double p_star_R);
}
#endif