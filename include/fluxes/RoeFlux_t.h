#ifndef FVM_ROEFLUX_T_H
#define FVM_ROEFLUX_T_H

#include "entities/FluxCalculator_t.h"
#include <vector>

namespace FVM { namespace Fluxes {
    class RoeFlux_t final : public FVM::Entities::FluxCalculator_t { 
        public: 
            RoeFlux_t(int n_faces);
            virtual ~RoeFlux_t();

            void invert_matrix(const double (&input)[9], double (&output)[9]);
            void multiply_matrix(const double (&left)[9], const double (&right)[9], double (&result)[9]);
            virtual void calculate_fluxes(double delta_t, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &F_1, std::vector<double> &F_2, std::vector<double> &F_3);
    };
}}
#endif