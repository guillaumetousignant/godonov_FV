#ifndef HLLEFLUX_T_H
#define HLLEFLUX_T_H

#include "FluxCalculator_t.h"
#include "Mesh1D_t.h"
#include <vector>

class HLLEFlux_t final : public FluxCalculator_t { 
    public: 
        HLLEFlux_t(int n_faces);
        virtual ~HLLEFlux_t();

        virtual void calculate_fluxes(Mesh1D_t &mesh, double delta_t);
};

#endif