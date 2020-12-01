// MCG 5136 - Finite-Volume Methods Assignment 3
// Guillaume Tousignant, 0300151859
// November 27th, 2020

#include <chrono>
#include <vector>
#include <iostream>
#include "Problem_t.h"

int main(void) {
    constexpr int n_problems = 4;
    constexpr int n_resolutions = 4;
    constexpr int n_points_analytic = 1000;

    const int n_points[] = {100,
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
    
    const double time[] = {12.0e-3,
                           25.0e-3,
                           35.0e-3,
                           7.0e-3};

    const double discontinuity[] = {2.0,
                                    2.0,
                                    5.0,
                                    5.0};
    
    std::vector<Problem_t> problems;
    problems.reserve(4);

    for (int i = 0; i < n_problems; ++i) {
        problems.push_back(Problem_t(rho[i][0], rho[i][1], u[i][0], u[i][1], p[i][0], p[i][1], time[i], discontinuity[i], n_points_analytic, i + 1));
    }

    // Starting actual computation
    auto t_start = std::chrono::high_resolution_clock::now();
    for (auto &problem : problems) {
        problem.solve();
    }
    auto t_end = std::chrono::high_resolution_clock::now();

    std::cout << "Computation time: " 
            << std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0 
            << "s." << std::endl;

    for (auto &problem : problems) {
        problem.write_solution();
    }

    return 0;
}