// MCG 5136 - Finite-Volume Methods Assignment 6
// Guillaume Tousignant, 0300151859
// December 21th, 2020

#include "FVM.h"
#include <chrono>
#include <vector>
#include <iostream>
#include <filesystem>

using FVM::Entities::Vec2f;

int main(int argc, char *argv[]) {
    std::string filename;
    if (argc < 2){
        std::cout << "Enter mesh path:" << std::endl;
        while (filename.empty()){
            std::getline(std::cin, filename);
        }
        if (filename == "q"){
            return -1;
        }
    }
    else{
        filename = argv[1];
    }

    std::filesystem::path filepath(filename);

    constexpr int n_problems = 1;
    constexpr double cfl = 0.5;

    double t_end = 0.5;

    state initial_conditions {  
        1.225,
        Vec2f(200, 0.0),
        101.325e3,
        1.4
    };

    std::filesystem::path output_path = "data/problem_4.dat";

    int problem_number = 4;

    std::cout << "Creating mesh... ";
    auto timer_start_mesh = std::chrono::high_resolution_clock::now();
    FVM::Entities::Mesh2D_t mesh(filepath);
    auto timer_end_mesh = std::chrono::high_resolution_clock::now();
        std::cout << "Took " 
            << std::chrono::duration<double, std::milli>(timer_end_mesh - timer_start_mesh).count()/1000.0 
            << "s." << std::endl;

    FVM::Solvers::GodonovSolver2D_t<FVM::Fluxes::HLLEFlux_t, FVM::Limiters::VenkatakrishnanLimiter_t> solver;

    std::cout << "Processing problem #" << problem_number << std::endl;

    std::cout << '\t' << "Setting initial conditions... ";
    auto timer_start = std::chrono::high_resolution_clock::now();
    mesh.initial_conditions(initial_conditions);
    auto timer_end = std::chrono::high_resolution_clock::now();
    std::cout << "Took " 
        << std::chrono::duration<double, std::milli>(timer_end - timer_start).count()/1000.0 
        << "s." << std::endl;

    std::cout << '\t' << "Solving... ";
    timer_start = std::chrono::high_resolution_clock::now();
    solver.solve(t_end, cfl, mesh);
    timer_end = std::chrono::high_resolution_clock::now();
    std::cout << "Took " 
        << std::chrono::duration<double, std::milli>(timer_end - timer_start).count()/1000.0 
        << "s." << std::endl;

    std::cout << '\t' << "Writing output file... ";
    timer_start = std::chrono::high_resolution_clock::now();
    mesh.write_tecplot(output_path, problem_number, t_end);
    timer_end = std::chrono::high_resolution_clock::now();
    std::cout << "Took " 
        << std::chrono::duration<double, std::milli>(timer_end - timer_start).count()/1000.0 
        << "s." << std::endl;

    return 0;
}