#ifndef FLUXCALCULATOR_T_H
#define FLUXCALCULATOR_T_H

#include <vector>

class FluxCalculator_t { 
    public: 
        FluxCalculator_t() {};
        virtual ~FluxCalculator_t() {};

        virtual void calculate_fluxes(double delta_t, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &F_1, std::vector<double> &F_2, std::vector<double> &F_3);
};

#endif