#ifndef ROEFLUX_T_H
#define ROEFLUX_T_H

#include "FluxCalculator_t.h"
#include "Mesh1D_t.h"
#include <vector>

class RoeFlux_t final : public FluxCalculator_t { 
    public: 
        RoeFlux_t(int n_faces);
        virtual ~RoeFlux_t();

        virtual void calculate_fluxes(Mesh1D_t &mesh, double delta_t);
};

#endif