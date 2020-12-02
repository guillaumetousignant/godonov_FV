#ifndef GODONOVSOLVER_T_H
#define GODONOVSOLVER_T_H

#include "Solver_t.h"
#include "Mesh1D_t.h"
#include <vector>

class GodonovSolver_t final : public Solver_t { 
    public: 
        GodonovSolver_t(double rho_L, double rho_R, double u_L, double u_R, double p_L, double p_R, double time, double discontinuity, int n_cells, int problem_number, double cfl);
        virtual ~GodonovSolver_t();

        Mesh1D_t mesh_;
        double cfl_;

        virtual void solve();

    private:
        double* a_star_L_; // These are split to be easier to deal with when working on a GPU
        double* a_star_R_;
        double* u_star_;
        double* p_star_L_;
        double* p_star_R_;
        double* p_star_prime_L_;
        double* p_star_prime_R_;
        double* C_L_;
        double* C_R_;

        double u_star_initial_guess(double a_L, double a_R, double u_L, double u_R, double p_L, double p_R, double gamma_L, double gamma_R);
        double calculate_C(double a, double p, double gamma);
        void calculate_a_star(double a_L, double a_R, double u_L, double u_R, double p_L, double p_R, double gamma_L, double gamma_R, double u_star, double a_star_L, double a_star_R, double p_star_L, double p_star_R);
        void left_shock(double a_L, double u_L, double p_L, double gamma_L, double C_L, double u_star, double &p_star_L, double &p_star_prime_L);
        void right_shock(double a_R, double u_R, double p_R, double gamma_R, double C_R, double u_star, double &p_star_R, double &p_star_prime_R);
        void left_rarefaction(double a_L, double u_L, double p_L, double gamma_L, double u_star, double &a_star_L, double &p_star_L, double &p_star_prime_L);
        void right_rarefaction(double a_R, double u_R, double p_R, double gamma_R, double u_star, double &a_star_R, double &p_star_R, double &p_star_prime_R);
        void solve_flux(double a_L, double a_R, double u_L, double u_R, double p_L, double p_R, double gamma_L, double gamma_R, double C_L, double C_R, double &u_star, double a_star_L, double a_star_R, double &p_star_L, double &p_star_R, double &p_star_prime_L, double &p_star_prime_R);
        void calculate_fluxes();
};

#endif