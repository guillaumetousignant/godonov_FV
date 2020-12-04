#ifndef ROEFLUX_T_H
#define ROEFLUX_T_H

#include "FluxCalculator_t.h"
#include "Mesh1D_t.h"
#include <vector>

class RoeFlux_t final : public FluxCalculator_t { 
    public: 
        RoeFlux_t(int n_faces);
        virtual ~RoeFlux_t();

        void invert_matrix(const double (&input)[9], double (&output)[9]);
        void multiply_matrix(const double (&left)[9], const double (&right)[9], double (&result)[9]);
        virtual void calculate_fluxes(Mesh1D_t &mesh, double delta_t);
};

#endif