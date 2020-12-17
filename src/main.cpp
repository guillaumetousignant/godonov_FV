// MCG 5136 - Finite-Volume Methods Assignment 6
// Guillaume Tousignant, 0300151859
// December 21th, 2020

#include "FVM.h"
#include <chrono>
#include <vector>
#include <iostream>
#include <filesystem>

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

    FVM::Entities::Mesh2D_t mesh(std::filesystem::path(filename));

    return 0;
}