#ifndef EXACTRIEMANNFLUX_T_H
#define EXACTRIEMANNFLUX_T_H

#include "FluxCalculator_t.h"
#include "Mesh1D_t.h"

class ExactRiemannFlux_t final : public FluxCalculator_t { 
    public: 
        ExactRiemannFlux_t(int n_faces);
        virtual ~ExactRiemannFlux_t();

        double* a_star_L_; // These are split to be easier to deal with when working on a GPU
        double* a_star_R_;
        double* u_star_;
        double* p_star_L_;
        double* p_star_R_;
        double* p_star_prime_L_;
        double* p_star_prime_R_;
        double* C_L_;
        double* C_R_;

        virtual void calculate_fluxes(Mesh1D_t &mesh, double delta_t);
};

#endif