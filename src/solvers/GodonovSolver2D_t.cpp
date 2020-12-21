#include "solvers/GodonovSolver2D_t.h"
#include "fluxes/ExactRiemannFlux_t.h"
#include "fluxes/RoeFlux_t.h"
#include "fluxes/RoeEntropyFlux_t.h"
#include "fluxes/HLLEFlux_t.h"
#include <cmath>
#include <algorithm>

using FVM::Solvers::GodonovSolver2D_t;
using FVM::Fluxes::ExactRiemannFlux_t;
using FVM::Fluxes::RoeFlux_t;
using FVM::Fluxes::RoeEntropyFlux_t;
using FVM::Fluxes::HLLEFlux_t;

template class GodonovSolver2D_t<ExactRiemannFlux_t>; // Like, I understand why I need this, but man is it crap.
template class GodonovSolver2D_t<RoeFlux_t>;
template class GodonovSolver2D_t<RoeEntropyFlux_t>;
template class GodonovSolver2D_t<HLLEFlux_t>;

template<typename FluxCalculator>
GodonovSolver2D_t<FluxCalculator>::GodonovSolver2D_t() :
        flux_calculator_() {}

template<typename FluxCalculator>
GodonovSolver2D_t<FluxCalculator>::~GodonovSolver2D_t() {}

template<typename FluxCalculator>
void GodonovSolver2D_t<FluxCalculator>::solve(double time, double cfl, FVM::Entities::Mesh2D_t& mesh) {
    double time = 0.0;

    while (time < end_time_) {
        double delta_t = calculate_delta_t(cfl, mesh);
        if (time + delta_t > end_time_) {
            delta_t = end_time_ - time;
        }

        mesh.boundary_conditions();
        flux_calculator_.calculate_fluxes(delta_t, mesh);
        timestep(delta_t, mesh);
        time += delta_t;
    }
}

template<typename FluxCalculator>
double GodonovSolver2D_t<FluxCalculator>::calculate_delta_t(double cfl, FVM::Entities::Mesh2D_t& mesh) {
    double min_dt = std::numeric_limits<double>::infinity();
    for (int i = 0; i <= mesh.n_cells_; ++i) {
        min_dt = std::min(min_dt, cfl * std::sqrt(mesh.cells_[i].area_) / std::abs(std::sqrt(std::pow(mesh.cells_[i].u_.x(), 2), std::pow(mesh.cells_[i].u_.y(), 2)) + mesh_.a_[i]));
    }
    return min_dt;
}

template<typename FluxCalculator>
void GodonovSolver2D_t<FluxCalculator>::timestep(double delta_t, FVM::Entities::Mesh2D_t& mesh) {
    
}