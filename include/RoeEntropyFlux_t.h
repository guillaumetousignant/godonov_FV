#ifndef ROEENTROPYFLUX_T_H
#define ROEENTROPYFLUX_T_H

#include "FluxCalculator_t.h"
#include "Mesh1D_t.h"
#include <vector>

class RoeEntropyFlux_t final : public FluxCalculator_t { 
    public: 
        RoeEntropyFlux_t(int n_faces);
        virtual ~RoeEntropyFlux_t();

        virtual void calculate_fluxes(Mesh1D_t &mesh, double delta_t);
};

#endif