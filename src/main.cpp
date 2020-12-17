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

    std::vector<Vec2f> centers{
        {0, 0},
        {0.3, 0.3},
        {0, 0}
    };

    std::vector<std::vector<state>> initial_conditions{
        {   
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

    std::vector<FVM::Entities::Mesh2D_t> meshes{
        filepath,
        filepath,
        filepath
    };

    for (size_t i = 0; i < initial_conditions.size(); ++i) {
        meshes[i].initial_conditions(centers[i], initial_conditions[i][0], initial_conditions[i][1], initial_conditions[i][2], initial_conditions[i][3]);
    }


    return 0;
}