#ifndef FVM_EXACTRIEMANNFLUX_T_H
#define FVM_EXACTRIEMANNFLUX_T_H

#include "entities/FluxCalculator_t.h"
#include <vector>

namespace FVM { namespace Fluxes {
    class ExactRiemannFlux_t final : public FVM::Entities::FluxCalculator_t { 
        public: 
            ExactRiemannFlux_t();
            virtual ~ExactRiemannFlux_t();

            virtual void calculate_fluxes(double delta_t, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &F_1, std::vector<double> &F_2, std::vector<double> &F_3);
    };
}}
#endif