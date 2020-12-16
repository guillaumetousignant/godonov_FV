#include "entities/Mesh2D_t.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

FVM::Entities::Mesh2D_t::Mesh2D_t(int n_cells, double delta_x) :
        n_cells_(n_cells),
        a_(n_cells_ + 2),  // Plus 2, because the ending cells are ghost cells
        u_(n_cells_ + 2),
        p_(n_cells_ + 2),
        x_(n_cells_ + 2),
        gamma_(n_cells_ + 2),
        delta_x_(delta_x),
        n_faces_(n_cells + 1),
        F_1_(n_faces_),
        F_2_(n_faces_),
        F_3_(n_faces_) {}

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

void FVM::Entities::Mesh2D_t::readSU2(const std::string &filename){
    std::string line;
    std::string token;
    int value;

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

    n_points_ = value;
    points_ = std::vector<Vec2f>(n_points_);

    for (unsigned int i = 0; i < n_points_; ++i) {
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

    n_cells_ = value;





    

    n_normals_ = n_points_;
    elements_ = new unsigned int[3 * n_elements_];
    element_normals_ = new unsigned int[3 * n_elements_];
    for (unsigned int i = 0; i < n_elements_; ++i) {
        std::getline(meshfile, line);
        std::istringstream liness2(line);
        liness2 >> token;
        if (token != "5") {
            std::cerr << "Error: expected token '5', found '" << token << "'. Exiting." << std::endl;
            return;
        }
        unsigned int val0, val1, val2;
        liness2 >> val0 >> val1 >> val2;
        elements_[3 * i] = val0 - 1;
        elements_[3 * i + 1] = val1 - 1;
        elements_[3 * i + 2] = val2 - 1;
        element_normals_[3 * i] = elements_[3 * i];
        element_normals_[3 * i + 1] = elements_[3 * i + 1];
        element_normals_[3 * i + 2] = elements_[3 * i + 2];
    }

    unsigned int n_markers;
    std::getline(meshfile, line);
    std::istringstream liness4(line);
    liness4 >> token;
    if (token == "NMARK="){
        liness4 >> n_markers;
    }
    else {
        std::cerr << "Error: expected marker 'NMARK=', found '" << token << "'. Exiting." << std::endl;
        return;
    }

    n_walls_ = n_markers - 1;
    n_wall_ = new unsigned int[n_walls_];
    unsigned int wall_index = 0;
    walls_ = new unsigned int*[n_walls_];
    sum_n_wall_ = 0;

    for (unsigned int i = 0; i < n_markers; ++i) {
        std::string type;
        std::getline(meshfile, line);
        std::istringstream liness5(line);
        liness5 >> token;
        if (token == "MARKER_TAG="){
            liness5 >> type;
        }
        else {
            std::cerr << "Error: expected marker 'MARKER_TAG=', found '" << token << "'. Exiting." << std::endl;
            return;
        }

        if (type == "wall") {
            std::getline(meshfile, line);
            std::istringstream liness2(line);
            liness2 >> token;
            if (token == "MARKER_ELEMS="){
                liness2 >> n_wall_[wall_index];
            }
            else {
                std::cerr << "Error: expected marker 'MARKER_ELEMS=', found '" << token << "'. Exiting." << std::endl;
                return;
            }

            sum_n_wall_ += n_wall_[wall_index];
            walls_[wall_index] = new unsigned int[2 * n_wall_[wall_index]];
            for (unsigned int j = 0; j < n_wall_[wall_index]; ++j) {

                std::getline(meshfile, line);
                std::istringstream liness6(line);

                liness6 >> token;
                if (token != "3") {
                    std::cerr << "Error: expected token '3', found '" << token << "'. Exiting." << std::endl;
                    return;
                }

                unsigned int val0, val1;
                liness6 >> val0 >> val1;
                walls_[wall_index][2 * j] = val0 - 1;
                walls_[wall_index][2 * j + 1] = val1 - 1;
            }
            ++wall_index;
        }
        else if (type == "farfield") {
            std::getline(meshfile, line);
            std::istringstream liness2(line);
            liness2 >> token;
            if (token == "MARKER_ELEMS="){
                liness2 >> n_farfield_;
            }
            else {
                std::cerr << "Error: expected marker 'MARKER_ELEMS=', found '" << token << "'. Exiting." << std::endl;
                return;
            }

            farfield_ = new unsigned int[2 * n_farfield_];
            for (unsigned int j = 0; j < n_farfield_; ++j) {

                std::getline(meshfile, line);
                std::istringstream liness6(line);

                liness6 >> token;
                if (token != "3") {
                    std::cerr << "Error: expected token '3', found '" << token << "'. Exiting." << std::endl;
                    return;
                }

                unsigned int val0, val1;
                liness6 >> val0 >> val1;
                farfield_[2 * j] = val0 - 1;
                farfield_[2 * j + 1] = val1 - 1;
            }
        }
        else {
            std::cerr << "Error: expected marker tag 'wall' or 'farfield', found '" << type << "'. Exiting." << std::endl;
            return;
        }
    }
    meshfile.close();

    normals_ = new Vec3f[n_points_];

    computeNodeToFace();
    computeNormals(n_points_);
}