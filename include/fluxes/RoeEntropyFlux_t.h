#ifndef FVM_ROEENTROPYFLUX_T_H
#define FVM_ROEENTROPYFLUX_T_H

#include "entities/FluxCalculator_t.h"
#include <vector>

namespace FVM { namespace Fluxes {
    class RoeEntropyFlux_t final : public FVM::Entities::FluxCalculator_t { 
        public: 
            RoeEntropyFlux_t();
            virtual ~RoeEntropyFlux_t();

            void invert_matrix(const double (&input)[9], double (&output)[9]);
            void multiply_matrix(const double (&left)[9], const double (&right)[9], double (&result)[9]);
            virtual void calculate_fluxes(double delta_t, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &F_1, std::vector<double> &F_2, std::vector<double> &F_3) override;
            virtual void calculate_fluxes_higher_order(double delta_t, const std::vector<double> x, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &F_1, std::vector<double> &F_2, std::vector<double> &F_3, const std::vector<double> du_dx, const std::vector<double> da_dx, const std::vector<double> dp_dx) override;
            virtual void calculate_fluxes(double delta_t, FVM::Entities::Mesh2D_t &mesh) override;
            virtual void calculate_fluxes_higher_order(double delta_t, FVM::Entities::Mesh2D_t &mesh) override;
            virtual void calculate_fluxes_higher_order_hat(double delta_t, FVM::Entities::Mesh2D_t &mesh) override;
    };
}}
#endif