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

}

template<typename FluxCalculator>
double GodonovSolver2D_t<FluxCalculator>::calculate_delta_t() {

    return 0;
}

template<typename FluxCalculator>
void GodonovSolver2D_t<FluxCalculator>::timestep(double delta_t, double delta_x, const std::vector<double> &gamma, std::vector<double> &u, std::vector<double> &a, std::vector<double> &p, const std::vector<double> F_1, const std::vector<double> F_2, const std::vector<double> F_3) {

}