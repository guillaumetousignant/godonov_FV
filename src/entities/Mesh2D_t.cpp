#include "entities/Mesh2D_t.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <filesystem>

FVM::Entities::Mesh2D_t::Mesh2D_t(std::filesystem::path filename) {
    if (filename.extension().string() == ".su2") {
        readSU2(filename);
    }
    else {
        std::cerr << "Error: file '" << filename << "' has extension '" << filename.extension() << "'. Only su2 is supported for now. Exiting." << std::endl;
        return; 
    }
}

FVM::Entities::Mesh2D_t::~Mesh2D_t() {}

void FVM::Entities::Mesh2D_t::initial_conditions(double a_L, double a_R, double u_L, double u_R, double p_L, double p_R, double x_L, double x_R, double gamma_L, double gamma_R, double discontinuity) {
    for (int i = 0; i < n_cells_ + 2; ++i) { // No need for explicit treatment of the boundary conditions, nice
        x_[i] = x_L + (i - 0.5) * delta_x_;
        if (x_[i] < discontinuity) {
            a_[i] = a_L;
            u_[i] = u_L;
            p_[i] = p_L;
            gamma_[i] = gamma_L;
        }
        else {
            a_[i] = a_R;
            u_[i] = u_R;
            p_[i] = p_R;
            gamma_[i] = gamma_R;
        }
    }
}

void FVM::Entities::Mesh2D_t::readSU2(std::filesystem::path filename){
    std::string line;
    std::string token;
    size_t value;

    std::ifstream meshfile(filename);
    if (!meshfile.is_open()) {
        std::cerr << "Error: file '" << filename << "' could not be opened. Exiting." << std::endl;
        return;
    }

    do {
        std::getline(meshfile, line);
        if (~line.empty()) {
            std::istringstream liness(line);
            liness >> token;
            liness >> value;
        }   
    }
    while (token != "NPOIN=");

    points_ = std::vector<Vec2f>(value);

    for (size_t i = 0; i < points_.size(); ++i) {
        std::getline(meshfile, line);
        std::istringstream liness2(line);
        liness2 >> points_[i][0] >> points_[i][1];
    }

    do {
        std::getline(meshfile, line);
        if (~line.empty()) {
            std::istringstream liness(line);
            liness >> token;
            liness >> value;
        }   
    }
    while (token != "NELEM=");

    cells_ = std::vector<Element_t>(value);
    n_cells_ = value;

    for (size_t i = 0; i < cells_.size(); ++i) {
        int n_sides;
        size_t val[4];

        std::getline(meshfile, line);
        std::istringstream liness2(line);
        liness2 >> token;
        if (token == "9") {
            n_sides = 4;
            liness2 >> val[0] >> val[1] >> val[2] >> val[3];

            cells_[i] = Element_t(n_sides);
            for (int j = 0; j < n_sides; ++j) {
                cells_[i].nodes_[j] = val[j] - 1;
            }
        }
        else {
            std::cerr << "Error: expected token '9', found '" << token << "'. Exiting." << std::endl;
            return;
        }
    }

    int n_markers;
    do {
        std::getline(meshfile, line);
        if (~line.empty()) {
            std::istringstream liness(line);
            liness >> token;
            liness >> n_markers;
        }   
    }
    while (token != "NMARK=");

    

    if (n_markers != 1) {
        std::cerr << "Error: expected 1 boundary, found '" << n_markers << "'. Exiting." << std::endl;
        return;
    }

    std::string type;
    do {
        std::getline(meshfile, line);
        if (~line.empty()) {
            std::istringstream liness(line);
            liness >> token;
            liness >> type;
        }   
    }
    while (token != "MARKER_TAG=");

    if (type == "farfield") {
        do {
            std::getline(meshfile, line);
            if (~line.empty()) {
                std::istringstream liness(line);
                liness >> token;
                liness >> value;
            }   
        }
        while (token != "MARKER_ELEMS=");

        n_boundary_ = value;
        cells_.reserve(n_cells_ + n_boundary_);

        for (size_t j = 0; j < n_boundary_; ++j) {

            std::getline(meshfile, line);
            std::istringstream liness6(line);

            liness6 >> token;
            if (token != "3") {
                std::cerr << "Error: expected token '3', found '" << token << "'. Exiting." << std::endl;
                return;
            }

            size_t val0, val1;
            liness6 >> val0 >> val1;
            cells_.push_back(Element_t(2));
            cells_[n_cells_ + j].nodes_[0] = val0 - 1;
            cells_[n_cells_ + j].nodes_[1] = val1 - 1;
        }
    }
    else {
        std::cerr << "Error: expected marker tag 'farfield', found '" << type << "'. Exiting." << std::endl;
        return;
    }

    meshfile.close();
}