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

    constexpr int n_problems = 3;

    double t_end[3] {
        0.75e-3,
        2.53e-3,
        1.9e-3
    };

    Vec2f centers[3] {
        {0, 0},
        {0.3, 0.3},
        {0, 0}
    };

    state initial_conditions[3][4] {
        {   
            state{
                1.225,
                Vec2f(0.0, 0.0),
                101.325e3,
                1.4
            },
            state{
                1.225,
                Vec2f(0.0, 0.0),
                101.325e3,
                1.4
            },
            state{
                0.30625,
                Vec2f(0.0, 0.0),
                25.33125e3,
                1.4
            },
            state{
                1.225,
                Vec2f(0.0, 0.0),
                101.325e3,
                1.4
            }
        },
        {   
            state{
                1.5,
                Vec2f(0.0, 0.0),
                150e3,
                1.4
            },
            state{
                0.532258064516129,
                Vec2f(381.385035698, 0.0),
                30e3,
                1.4
            },
            state{
                0.137992831541219,
                Vec2f(381.385035698, 381.385035698),
                2.9032258064516e3,
                1.4
            },
            state{
                0.532258064516129,
                Vec2f(0.0,  381.385035698),
                30e3,
                1.4
            }
        },
        {   
            state{
                1,
                Vec2f(237.170824513, -158.113883008),
                100e3,
                1.4
            },
            state{
                3,
                Vec2f(237.170824513, 158.113883008),
                100e3,
                1.4
            },
            state{
                1,
                Vec2f(-237.170824513, 158.113883008),
                100e3,
                1.4
            },
            state{
                3,
                Vec2f(-237.170824513, -158.113883008),
                100e3,
                1.4
            }
        },
    };

    std::filesystem::path output_paths[3] {
        "data/problem_1.dat",
        "data/problem_2.dat",
        "data/problem_3.dat"
    };

    int problem_numbers[3] {
        1,
        2,
        3
    };

    FVM::Entities::Mesh2D_t meshes[3] {
        filepath,
        filepath,
        filepath
    };

    for (int i = 0; i < n_problems; ++i) {
        meshes[i].initial_conditions(centers[i], initial_conditions[i][0], initial_conditions[i][1], initial_conditions[i][2], initial_conditions[i][3]);
    }

    for (int i = 0; i < n_problems; ++i) {
        meshes[i].write_tecplot(output_paths[i], problem_numbers[i], t_end[i]);
    };

    return 0;
}