#ifndef ROEENTROPYFLUX_T_H
#define ROEENTROPYFLUX_T_H

#include "FluxCalculator_t.h"
#include "Mesh1D_t.h"
#include <vector>

class RoeEntropyFlux_t final : public FluxCalculator_t { 
    public: 
        RoeEntropyFlux_t(int n_faces);
        virtual ~RoeEntropyFlux_t();

        void invert_matrix(const double (&input)[9], double (&output)[9]);
        void multiply_matrix(const double (&left)[9], const double (&right)[9], double (&result)[9]);
        virtual void calculate_fluxes(Mesh1D_t &mesh, double delta_t);
};

#endif