#ifndef FVM_FLUXCALCULATOR_T_H
#define FVM_FLUXCALCULATOR_T_H

#include "entities/Mesh2D_t.h"
#include <vector>

namespace FVM { namespace Entities {
    class FluxCalculator_t { 
        public: 
            FluxCalculator_t() {};
            virtual ~FluxCalculator_t() {};

            virtual void calculate_fluxes(double delta_t, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &F_1, std::vector<double> &F_2, std::vector<double> &F_3) = 0;
            virtual void calculate_fluxes_higher_order(double delta_t, const std::vector<double> x, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &F_1, std::vector<double> &F_2, std::vector<double> &F_3, const std::vector<double> du_dx, const std::vector<double> da_dx, const std::vector<double> dp_dx) = 0;
            virtual void calculate_fluxes(double delta_t, FVM::Entities::Mesh2D_t &mesh) = 0;
            virtual void calculate_fluxes_higher_order(double delta_t, FVM::Entities::Mesh2D_t &mesh) = 0;
    };
}}
#endif