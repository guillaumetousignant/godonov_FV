#include "entities/Mesh2D_t.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <filesystem>
#include <cmath>

using FVM::Entities::Vec2f;

constexpr double R = 8.31446261815324;

FVM::Entities::Mesh2D_t::Mesh2D_t(std::filesystem::path filename) {
    if (filename.extension().string() == ".su2") {
        read_su2(filename);
        compute_cell_geometry();
        build_node_to_cell();
        build_cell_to_cell();
        build_faces();
        compute_face_geometry();
    }
    else {
        std::cerr << "Error: file '" << filename << "' has extension '" << filename.extension() << "'. Only su2 is supported for now. Exiting." << std::endl;
        return; 
    }
}

FVM::Entities::Mesh2D_t::~Mesh2D_t() {}

void FVM::Entities::Mesh2D_t::initial_conditions(FVM::Entities::Vec2f center, const state& state_NE, const state& state_NW, const state& state_SW, const state& state_SE) {
    for (auto& cell: cells_) {
        if (cell.center_.x() > center.x() && cell.center_.y() >= center.y()) {
            cell.a_ = std::sqrt(state_NE.gamma * state_NE.p / state_NE.rho);
            cell.u_ = state_NE.u;
            cell.p_ = state_NE.p;
            cell.gamma_ = state_NE.gamma;
        }
        else if (cell.center_.x() <= center.x() && cell.center_.y() > center.y()) {
            cell.a_ = std::sqrt(state_NW.gamma * state_NW.p / state_NW.rho);
            cell.u_ = state_NW.u;
            cell.p_ = state_NW.p;
            cell.gamma_ = state_NW.gamma;
        }
        else if (cell.center_.x() < center.x() && cell.center_.y() <= center.y()) {
            cell.a_ = std::sqrt(state_SW.gamma * state_SW.p / state_SW.rho);
            cell.u_ = state_SW.u;
            cell.p_ = state_SW.p;
            cell.gamma_ = state_SW.gamma;
        }
        else {
            cell.a_ = std::sqrt(state_SE.gamma * state_SE.p / state_SE.rho);
            cell.u_ = state_SE.u;
            cell.p_ = state_SE.p;
            cell.gamma_ = state_SE.gamma;
        }
    }
}

void FVM::Entities::Mesh2D_t::boundary_conditions() {
    for (size_t i = n_cells_; i < n_cells_ + n_boundary_; ++i) {
        // CHECK this will not work with 2nd order
        cells_[i].a_ = cells_[cells_[i].cells_[0]].a_; // Both cell to cell links point to the cell inside the domain for boundary line elements.
        cells_[i].u_ = cells_[cells_[i].cells_[0]].u_;
        cells_[i].p_ = cells_[cells_[i].cells_[0]].p_;
        cells_[i].gamma_ = cells_[cells_[i].cells_[0]].gamma_;
    }
}

void FVM::Entities::Mesh2D_t::write_tecplot(std::filesystem::path filename, int problem_number, double time) {
    std::ofstream file;

    std::filesystem::create_directory(filename.parent_path());

    file.open(filename);

    file << "TITLE = \"Problem " << problem_number << " at t= " << time << "\"" << std::endl;
    file << "VARIABLES = \"X\", \"Y\", \"U_x\", \"U_y\", \"rho\", \"p\", \"mach\", \"T\"" << std::endl;
    file << "ZONE T=\"Zone 1\", ZONETYPE=FEQUADRILATERAL, NODES=" << nodes_.size() << ", ELEMENTS=" << n_cells_ << ", DATAPACKING=BLOCK, VARLOCATION=([1,2]=nodal,[3,4,5,6,7,8]=cellcentered), SOLUTIONTIME=" << time << std::endl;

    for (const auto& node: nodes_) {
        file << std::setw(12) << node.pos_.x() << " ";
    }
    file << std::endl;

    for (const auto& node: nodes_) {
        file << std::setw(12) << node.pos_.y() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].u_.x() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].u_.y() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].gamma_ * cells_[i].p_/std::pow(cells_[i].a_, 2) << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].p_ << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << std::sqrt(std::pow(cells_[i].u_.x(), 2) + std::pow(cells_[i].u_.y(), 2)) / std::sqrt(cells_[i].p_ / std::pow(cells_[i].gamma_ * cells_[i].p_/std::pow(cells_[i].a_, 2), cells_[i].gamma_)) << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].p_/(cells_[i].gamma_ * cells_[i].p_/std::pow(cells_[i].a_, 2) * R) << " ";
    }
    file << std::endl;

    // Connectivity
    file << std::endl;
    for (int i = 0; i < n_cells_; ++i) {
        if (cells_[i].nodes_.size() < 4) {
            for (int j = 0; j < 4 - cells_[i].nodes_.size(); ++j) {
                file << cells_[i].nodes_[0] + 1 << " ";
            }
        }
        for (int j = 0; j < cells_[i].nodes_.size(); ++j) {
            file << cells_[i].nodes_[j] + 1 << " ";
        }
        file << std::endl;
    }

    file.close();
}

void FVM::Entities::Mesh2D_t::read_su2(std::filesystem::path filename){
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
        if (!line.empty()) {
            std::istringstream liness(line);
            liness >> token;
            liness >> value;
        }   
    }
    while (token != "NPOIN=");

    nodes_ = std::vector<Node_t>(value);

    for (size_t i = 0; i < nodes_.size(); ++i) {
        std::getline(meshfile, line);
        std::istringstream liness2(line);
        liness2 >> nodes_[i].pos_[0] >> nodes_[i].pos_[1];
    }

    do {
        std::getline(meshfile, line);
        if (!line.empty()) {
            std::istringstream liness(line);
            liness >> token;
            liness >> value;
        }   
    }
    while (token != "NELEM=");

    cells_ = std::vector<Cell_t>(value);
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

            cells_[i] = Cell_t(n_sides);
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
        if (!line.empty()) {
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
        if (!line.empty()) {
            std::istringstream liness(line);
            liness >> token;
            liness >> type;
        }   
    }
    while (token != "MARKER_TAG=");

    if (type == "farfield") {
        do {
            std::getline(meshfile, line);
            if (!line.empty()) {
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
            cells_.push_back(Cell_t(2));
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

void FVM::Entities::Mesh2D_t::build_node_to_cell() {
    for (size_t i = 0; i < nodes_.size(); ++i) {
        for (size_t j = 0; j < cells_.size(); ++j) {
            for (int k = 0; k < cells_[j].nodes_.size(); ++k) {
                if (cells_[j].nodes_[k] == i) {
                    nodes_[i].cells_.push_back(j);
                    break;
                }
            }
        }
    }
}

void FVM::Entities::Mesh2D_t::build_cell_to_cell() {
    for (size_t i = 0; i < cells_.size(); ++i) {
        for (size_t j = 0; j < cells_[i].nodes_.size() - 1; ++j) {
            for (size_t m = 0; m < nodes_[cells_[i].nodes_[j]].cells_.size(); ++m) {
                if (nodes_[cells_[i].nodes_[j]].cells_[m] != i) {
                    for (size_t n = 0; n < nodes_[cells_[i].nodes_[j + 1]].cells_.size(); ++n) {
                        if (nodes_[cells_[i].nodes_[j]].cells_[m] == nodes_[cells_[i].nodes_[j + 1]].cells_[n]) {
                            cells_[i].cells_[j] = nodes_[cells_[i].nodes_[j]].cells_[m];
                            goto endloop; // I hate this too don't worry
                        }
                    }
                }
            }
            endloop: ;
        }
        for (size_t m = 0; m < nodes_[cells_[i].nodes_[cells_[i].nodes_.size() - 1]].cells_.size(); ++m) {
            if (nodes_[cells_[i].nodes_.size() - 1].cells_[m] != i) {
                for (size_t n = 0; n < nodes_[cells_[i].nodes_[0]].cells_.size(); ++n) {
                    if (nodes_[cells_[i].nodes_[cells_[i].nodes_.size() - 1]].cells_[m] == nodes_[cells_[i].nodes_[0]].cells_[n]) {
                        cells_[i].cells_[cells_[i].nodes_.size() - 1] = nodes_[cells_[i].nodes_[cells_[i].nodes_.size() - 1]].cells_[m];
                        goto endloop2; // I hate this too don't worry
                    }
                }
            }
        }
        endloop2: ;
    }
}

void FVM::Entities::Mesh2D_t::build_faces() {
    size_t total_edges = 0;
    for (const auto& cell: cells_) {
        total_edges += cell.nodes_.size();
    }

    faces_.reserve(total_edges/2); // This is not exact

    for (size_t i = 0; i < cells_.size(); ++i) {
        for (size_t j = 0; j < cells_[i].nodes_.size() - 1; ++j) {
            size_t nodes[] = {cells_[i].nodes_[j], cells_[i].nodes_[j + 1]};
            bool found = false;
            for (size_t k = 0; k < faces_.size(); ++k) {
                if ((faces_[k].nodes_[0] == nodes[0] && faces_[k].nodes_[1] == nodes[1]) || (faces_[k].nodes_[0] == nodes[1] && faces_[k].nodes_[1] == nodes[0])) {
                    found = true;
                    faces_[k].cells_[1] = i;
                    cells_[i].faces_[j] = k;
                    break;
                }

            }

            if (!found) {
                cells_[i].faces_[j] = faces_.size();
                faces_.push_back(Face_t(nodes[0], nodes[1], i, -1));
            }
        }
        size_t nodes[] = {cells_[i].nodes_[cells_[i].nodes_.size() - 1], cells_[i].nodes_[0]};
        bool found = false;
        for (size_t k = 0; k < faces_.size(); ++k) {
            if ((faces_[k].nodes_[0] == nodes[0] && faces_[k].nodes_[1] == nodes[1]) || (faces_[k].nodes_[0] == nodes[1] && faces_[k].nodes_[1] == nodes[0])) {
                found = true;
                faces_[k].cells_[1] = i;
                cells_[i].faces_[cells_[i].nodes_.size() - 1] = k;
                break;
            }

        }

        if (!found) {
            cells_[i].faces_[cells_[i].nodes_.size() - 1] = faces_.size();
            faces_.push_back(Face_t(nodes[0], nodes[1], i, -1));
        }
    }
}

void FVM::Entities::Mesh2D_t::compute_cell_geometry() {
    for (auto& cell: cells_) {
        cell.center_ = Vec2f();
        for (auto node: cell.nodes_) {
            cell.center_ += nodes_[node].pos_;
        }
        cell.center_ /= cell.nodes_.size();

        // Area is calculated by triangularisation, concave shapes will generally NOT work.
        cell.area_ = 0; 
        for (size_t i = 0; i < cell.nodes_.size() - 2; ++i) {
            const Vec2f points[] = {nodes_[cell.nodes_[0]].pos_,
                                    nodes_[cell.nodes_[i + 1]].pos_,
                                    nodes_[cell.nodes_[i + 2]].pos_};
            cell.area_ += std::abs(points[0].x() * (points[1].y() - points[2].y()) + points[1].x() * (points[2].y() - points[0].y()) + points[2].x() * (points[0].y() - points[1].y())) * 0.5;
        }
    }
}

void FVM::Entities::Mesh2D_t::compute_face_geometry() {
    for (auto& face: faces_) {
        face.tangent_ = (nodes_[face.nodes_[1]].pos_ - nodes_[face.nodes_[0]].pos_).normalize(); // CHECK should be normalized or not?
        face.normal_ = Vec2f(face.tangent_.y(), -face.tangent_.x()); 

        const Vec2f center = (nodes_[face.nodes_[0]].pos_ + nodes_[face.nodes_[1]].pos_) * 0.5;
        const Vec2f delta = center - cells_[face.cells_[0]].center_;
        const double sign = std::copysign(1.0, face.normal_.dot(delta));
        face.normal_ *= sign;
        face.tangent_ *= sign;
    }
}
