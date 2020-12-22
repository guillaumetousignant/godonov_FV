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
void GodonovSolver2D_t<FluxCalculator>::solve(double end_time, double cfl, FVM::Entities::Mesh2D_t& mesh) {
    double time = 0.0;

    while (time < end_time) {
        double delta_t = calculate_delta_t(cfl, mesh);
        if (time + delta_t > end_time) {
            delta_t = end_time - time;
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
    for (size_t i = 0; i < mesh.n_cells_; ++i) {
        min_dt = std::min(min_dt, cfl * std::sqrt(mesh.cells_[i].area_) / std::abs(mesh.cells_[i].u_.magnitude() + mesh.cells_[i].a_));
    }
    return min_dt;
}

template<typename FluxCalculator>
void GodonovSolver2D_t<FluxCalculator>::timestep(double delta_t, FVM::Entities::Mesh2D_t& mesh) {
    //#pragma omp parallel for schedule(guided)
    for (long long i = 0; i < mesh.n_cells_; ++i) { // Because msvc openmp wants a signed index, so no size_t :(
        FVM::Entities::Cell_t& cell =  mesh.cells_[i];

        double F_1 = 0;
        double F_2 = 0;
        double F_3 = 0;
        double F_4 = 0;
        for (auto face_index: cell.faces_) {
            const FVM::Entities::Face_t& face = mesh.faces_[face_index];
            const FVM::Entities::Vec2f delta = face.center_ - cell.center_;
            const double sign = std::copysign(1.0, delta.dot(face.normal_));
            F_1 += face.F_1_.dot(face.normal_) * face.length_ * sign; // CHECK do this when computing fluxes, it's computed twice as it is now.
            F_2 += face.F_2_.dot(face.normal_) * face.length_ * sign;
            F_3 += face.F_3_.dot(face.normal_) * face.length_ * sign;
            F_4 += face.F_4_.dot(face.normal_) * face.length_ * sign;
        }

        const double U_1 = cell.gamma_ * cell.p_ /std::pow(cell.a_, 2) - delta_t * F_1/cell.area_;
        const double U_2 = cell.gamma_ * cell.p_ * cell.u_.x()/std::pow(cell.a_, 2) - delta_t * F_2/cell.area_;
        const double U_3 = cell.gamma_ * cell.p_ * cell.u_.y()/std::pow(cell.a_, 2) - delta_t * F_3/cell.area_;
        const double U_4 = cell.p_/(cell.gamma_ - 1) + cell.gamma_ * cell.p_ * cell.u_.magnitudeSquared() * 0.5/std::pow(cell.a_, 2) - delta_t * F_4/cell.area_;

        cell.u_ = FVM::Entities::Vec2f(U_2/U_1, U_3/U_1);
        cell.p_ = (cell.gamma_ - 1) * (U_4 - 0.5 * (U_2 * cell.u_.x() + U_3 * cell.u_.y()));
        cell.a_ = std::sqrt(cell.gamma_ * cell.p_ /U_1);
    }
}