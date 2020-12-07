#ifndef FVM_EXACTRIEMANNFLUX_T_H
#define FVM_EXACTRIEMANNFLUX_T_H

#include "entities/FluxCalculator_t.h"
#include <vector>

class ExactRiemannFlux_t final : public FluxCalculator_t { 
    public: 
        ExactRiemannFlux_t(int n_faces);
        virtual ~ExactRiemannFlux_t();

        std::vector<double> a_star_L_; // These are split to be easier to deal with when working on a GPU
        std::vector<double> a_star_R_;
        std::vector<double> u_star_;
        std::vector<double> p_star_L_;
        std::vector<double> p_star_R_;
        std::vector<double> p_star_prime_L_;
        std::vector<double> p_star_prime_R_;
        std::vector<double> C_L_;
        std::vector<double> C_R_;

        virtual void calculate_fluxes(double delta_t, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &F_1, std::vector<double> &F_2, std::vector<double> &F_3);
};

#endif