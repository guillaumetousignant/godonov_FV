// MCG 5136 - Finite-Volume Methods Assignment 3
// Guillaume Tousignant, 0300151859
// November 27th, 2020

#include "FVM.h"
#include <chrono>
#include <vector>
#include <iostream>

using FVM::Solvers::ExactSolver_t;
using FVM::Solvers::GodonovSolver_t;
using FVM::Solvers::GodonovSolverHigherOrder_t;
using FVM::Fluxes::ExactRiemannFlux_t;
using FVM::Fluxes::RoeFlux_t;
using FVM::Fluxes::RoeEntropyFlux_t;
using FVM::Fluxes::HLLEFlux_t;

int main(void) {
    constexpr int n_problems = 4;
    constexpr int n_resolutions = 4;
    constexpr int n_points_exact = 10000; // It's one "cell" so we need many points
    constexpr int n_points = 10; // 10 points per cell should be enough even for higher order
    constexpr double cfl = 0.5;

    const int n_cells[] = {100,
                            200,
                            500,
                            1000};

    const double rho[4][2] = {{2.281, 1.408},
                              {1.045, 3.483},
                              {1.598, 2.787},
                              {4.696, 1.408}};

    const double u[4][2] = {{164.83, 0.0},
                            {200, 200},
                            {-383.64, -216.97},
                            {0.0, 0.0}};

    const double p[4][2] = {{201.170e3, 101.1e3},
                            {300.0e3, 300.0e3},
                            {91.88e3, 200.0e3},
                            {404.4e3, 101.1e3}};
    
    const double end_time[] = {12.0e-3,
                               25.0e-3,
                               35.0e-3,
                               7.0e-3};

    const double discontinuity[] = {2.0,
                                    2.0,
                                    5.0,
                                    5.0};

    const double x_span[4][2] = {{0.0, 10.0},
                                 {0.0, 10.0},
                                 {0.0, 10.0},
                                 {0.0, 10.0}};
    
    std::vector<ExactSolver_t> exact_solution;
    std::vector<GodonovSolverHigherOrder_t<ExactRiemannFlux_t>> riemann_solution;
    std::vector<GodonovSolverHigherOrder_t<RoeFlux_t>> roe_solution;
    std::vector<GodonovSolverHigherOrder_t<RoeEntropyFlux_t>> roe_entropy_solution;
    std::vector<GodonovSolverHigherOrder_t<HLLEFlux_t>> hlle_solution;
    exact_solution.reserve(n_problems);
    riemann_solution.reserve(n_problems * n_resolutions);
    roe_solution.reserve(n_problems * n_resolutions);
    roe_entropy_solution.reserve(n_problems * n_resolutions);
    hlle_solution.reserve(n_problems * n_resolutions);

    for (int i = 0; i < n_problems; ++i) {
        exact_solution.push_back(ExactSolver_t(rho[i][0], rho[i][1], u[i][0], u[i][1], p[i][0], p[i][1], x_span[i][0], x_span[i][1], end_time[i], discontinuity[i], n_points_exact, i + 1));
        for (int j = 0; j < n_resolutions; ++j) {
            riemann_solution.push_back(GodonovSolverHigherOrder_t<ExactRiemannFlux_t>(rho[i][0], rho[i][1], u[i][0], u[i][1], p[i][0], p[i][1], x_span[i][0], x_span[i][1], end_time[i], discontinuity[i], n_points, n_cells[j], i + 1, cfl));
            roe_solution.push_back(GodonovSolverHigherOrder_t<RoeFlux_t>(rho[i][0], rho[i][1], u[i][0], u[i][1], p[i][0], p[i][1], x_span[i][0], x_span[i][1], end_time[i], discontinuity[i], n_points, n_cells[j], i + 1, cfl));
            roe_entropy_solution.push_back(GodonovSolverHigherOrder_t<RoeEntropyFlux_t>(rho[i][0], rho[i][1], u[i][0], u[i][1], p[i][0], p[i][1], x_span[i][0], x_span[i][1], end_time[i], discontinuity[i], n_points, n_cells[j], i + 1, cfl));
            hlle_solution.push_back(GodonovSolverHigherOrder_t<HLLEFlux_t>(rho[i][0], rho[i][1], u[i][0], u[i][1], p[i][0], p[i][1], x_span[i][0], x_span[i][1], end_time[i], discontinuity[i], n_points, n_cells[j], i + 1, cfl));
        }
    }

    // Starting computation
    auto t_start_exact = std::chrono::high_resolution_clock::now();
    for (auto &problem : exact_solution) {
        problem.solve();
    }
    auto t_end_exact = std::chrono::high_resolution_clock::now();

    std::cout << "Exact solution computation time: " 
            << std::chrono::duration<double, std::milli>(t_end_exact - t_start_exact).count()/1000.0 
            << "s." << std::endl;

    auto t_start_riemann = std::chrono::high_resolution_clock::now();
    for (auto &problem : riemann_solution) {
        problem.solve();
    }
    auto t_end_riemann = std::chrono::high_resolution_clock::now();

    std::cout << "Exact Riemann solver computation time: " 
            << std::chrono::duration<double, std::milli>(t_end_riemann - t_start_riemann).count()/1000.0 
            << "s." << std::endl;

    auto t_start_roe = std::chrono::high_resolution_clock::now();
    for (auto &problem : roe_solution) {
        problem.solve();
    }
    auto t_end_roe = std::chrono::high_resolution_clock::now();

    std::cout << "Roe solver computation time: " 
            << std::chrono::duration<double, std::milli>(t_end_roe - t_start_roe).count()/1000.0 
            << "s." << std::endl;

    auto t_start_roe_entropy = std::chrono::high_resolution_clock::now();
    for (auto &problem : roe_entropy_solution) {
        problem.solve();
    }
    auto t_end_roe_entropy = std::chrono::high_resolution_clock::now();

    std::cout << "Roe entropy fix solver computation time: " 
            << std::chrono::duration<double, std::milli>(t_end_roe_entropy - t_start_roe_entropy).count()/1000.0 
            << "s." << std::endl;

    auto t_start_hlle = std::chrono::high_resolution_clock::now();
    for (auto &problem : hlle_solution) {
        problem.solve();
    }
    auto t_end_hlle = std::chrono::high_resolution_clock::now();

    std::cout << "HLLE solver computation time: " 
            << std::chrono::duration<double, std::milli>(t_end_hlle - t_start_hlle).count()/1000.0 
            << "s." << std::endl;

    // Output
    for (auto &problem : exact_solution) {
        problem.write_solution("_exact");
    }

    for (auto &problem : riemann_solution) {
        problem.write_solution("_riemann_N" + std::to_string(problem.mesh_.n_cells_));
    }

    for (auto &problem : roe_solution) {
        problem.write_solution("_roe_N" + std::to_string(problem.mesh_.n_cells_));
    }
    for (auto &problem : roe_entropy_solution) {
        problem.write_solution("_roe_entropy_N" + std::to_string(problem.mesh_.n_cells_));
    }
    
    for (auto &problem : hlle_solution) {
        problem.write_solution("_hlle_N" + std::to_string(problem.mesh_.n_cells_));
    }

    return 0;
}