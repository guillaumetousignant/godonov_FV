#ifndef HLLEFLUX_T_H
#define HLLEFLUX_T_H

#include "FluxCalculator_t.h"
#include <vector>

class HLLEFlux_t final : public FluxCalculator_t { 
    public: 
        HLLEFlux_t(int n_faces);
        virtual ~HLLEFlux_t();

        virtual void calculate_fluxes(double delta_t, const std::vector<double> &gamma, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &F_1, std::vector<double> &F_2, std::vector<double> &F_3);
};

#endif