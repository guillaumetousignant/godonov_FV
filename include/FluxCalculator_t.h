#ifndef FLUXCALCULATOR_T_H
#define FLUXCALCULATOR_T_H

#include "Mesh1D_t.h"

class FluxCalculator_t { 
    public: 
        FluxCalculator_t() {};
        virtual ~FluxCalculator_t() {};

        virtual void calculate_fluxes(Mesh1D_t &mesh, double delta_t) = 0;
};

#endif