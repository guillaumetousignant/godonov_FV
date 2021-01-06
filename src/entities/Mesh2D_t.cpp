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
        exit(5); 
    }
}

void FVM::Entities::Mesh2D_t::initial_conditions(FVM::Entities::Vec2f center, const state& state_NE, const state& state_NW, const state& state_SW, const state& state_SE) {
    #pragma omp parallel for schedule(static)
    for (long long i = 0; i < cells_.size(); ++i) {
        FVM::Entities::Cell_t &cell = cells_[i];
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
    #pragma omp parallel for schedule(static)
    for (long long i = n_cells_; i < n_cells_ + n_farfield_; ++i) {
        // CHECK this will not work with 2nd order
        cells_[i].a_ = cells_[cells_[i].cells_[0]].a_; // Both cell to cell links point to the cell inside the domain for boundary line elements.
        cells_[i].u_ = cells_[cells_[i].cells_[0]].u_;
        cells_[i].p_ = cells_[cells_[i].cells_[0]].p_;
        cells_[i].gamma_ = cells_[cells_[i].cells_[0]].gamma_;
        cells_[i].a_derivative_ = {0, 0};
        cells_[i].ux_derivative_ = {0, 0};
        cells_[i].uy_derivative_ = {0, 0};
        cells_[i].p_derivative_ = {0, 0};
        cells_[i].gamma_derivative_ = {0, 0};
    }

    #pragma omp parallel for schedule(static)
    for (long long i = n_cells_ + n_farfield_; i < n_cells_ + n_farfield_ + n_wall_; ++i) {
        // CHECK this doesn't work with higher-order, as the velocity used to compute the flux are interpolated to the boundary using the derivative
        Vec2f u_prime(cells_[cells_[i].cells_[0]].u_.dot(faces_[cells_[i].faces_[0]].normal_), cells_[cells_[i].cells_[0]].u_.dot(faces_[cells_[i].faces_[0]].tangent_));
        u_prime.x() = -u_prime.x();

        const Vec2f normal_inv = Vec2f(faces_[cells_[i].faces_[0]].normal_.x(), faces_[cells_[i].faces_[0]].tangent_.x());
        const Vec2f tangent_inv = Vec2f(faces_[cells_[i].faces_[0]].normal_.y(), faces_[cells_[i].faces_[0]].tangent_.y());

        // CHECK this will not work with 2nd order
        cells_[i].a_ = cells_[cells_[i].cells_[0]].a_; // Both cell to cell links point to the cell inside the domain for boundary line elements.
        cells_[i].u_ = {normal_inv.dot(u_prime), tangent_inv.dot(u_prime)};
        cells_[i].p_ = cells_[cells_[i].cells_[0]].p_;
        cells_[i].gamma_ = cells_[cells_[i].cells_[0]].gamma_;
        cells_[i].a_derivative_ = {0, 0};
        cells_[i].ux_derivative_ = {0, 0};
        cells_[i].uy_derivative_ = {0, 0};
        cells_[i].p_derivative_ = {0, 0};
        cells_[i].gamma_derivative_ = {0, 0};
    }
}

void FVM::Entities::Mesh2D_t::boundary_conditions_hat() {
    #pragma omp parallel for schedule(static)
    for (long long i = n_cells_; i < n_cells_ + n_farfield_; ++i) {
        // CHECK this will not work with 2nd order
        cells_[i].a_hat_ = cells_[cells_[i].cells_[0]].a_hat_; // Both cell to cell links point to the cell inside the domain for boundary line elements.
        cells_[i].u_hat_ = cells_[cells_[i].cells_[0]].u_hat_;
        cells_[i].p_hat_ = cells_[cells_[i].cells_[0]].p_hat_;
        cells_[i].gamma_hat_ = cells_[cells_[i].cells_[0]].gamma_hat_;
        cells_[i].a_derivative_hat_ = {0, 0};
        cells_[i].ux_derivative_hat_ = {0, 0};
        cells_[i].uy_derivative_hat_ = {0, 0};
        cells_[i].p_derivative_hat_ = {0, 0};
        cells_[i].gamma_derivative_hat_ = {0, 0};
    }

    #pragma omp parallel for schedule(static)
    for (long long i = n_cells_ + n_farfield_; i < n_cells_ + n_farfield_ + n_wall_; ++i) {
        // CHECK this doesn't work with higher-order, as the velocity used to compute the flux are interpolated to the boundary using the derivative
        Vec2f u_prime(cells_[cells_[i].cells_[0]].u_hat_.dot(faces_[cells_[i].faces_[0]].normal_), cells_[cells_[i].cells_[0]].u_hat_.dot(faces_[cells_[i].faces_[0]].tangent_));
        u_prime.x() = -u_prime.x();

        const Vec2f normal_inv = Vec2f(faces_[cells_[i].faces_[0]].normal_.x(), faces_[cells_[i].faces_[0]].tangent_.x());
        const Vec2f tangent_inv = Vec2f(faces_[cells_[i].faces_[0]].normal_.y(), faces_[cells_[i].faces_[0]].tangent_.y());

        // CHECK this will not work with 2nd order
        cells_[i].a_hat_ = cells_[cells_[i].cells_[0]].a_hat_; // Both cell to cell links point to the cell inside the domain for boundary line elements.
        cells_[i].u_hat_ = {normal_inv.dot(u_prime), tangent_inv.dot(u_prime)};
        cells_[i].p_hat_ = cells_[cells_[i].cells_[0]].p_hat_;
        cells_[i].gamma_hat_ = cells_[cells_[i].cells_[0]].gamma_hat_;
        cells_[i].a_derivative_hat_ = {0, 0};
        cells_[i].ux_derivative_hat_ = {0, 0};
        cells_[i].uy_derivative_hat_ = {0, 0};
        cells_[i].p_derivative_hat_ = {0, 0};
        cells_[i].gamma_derivative_hat_ = {0, 0};
    }
}

void FVM::Entities::Mesh2D_t::write_tecplot(std::filesystem::path filename, int problem_number, double time) {
    std::ofstream file;

    std::filesystem::create_directory(filename.parent_path());

    file.open(filename);

    file << "TITLE = \"Problem " << problem_number << " at t= " << time << "\"" << std::endl;
    file << "VARIABLES = \"X\", \"Y\", \"U_x\", \"U_y\", \"rho\", \"p\", \"mach\", \"T\", \"da_dx\", \"da_dy\", \"dux_dx\", \"dux_dy\", \"duy_dx\", \"duy_dy\", \"dp_dx\", \"dp_dy\", \"da_dx_hat\", \"da_dy_hat\", \"dux_dx_hat\", \"dux_dy_hat\", \"duy_dx_hat\", \"duy_dy_hat\", \"dp_dx_hat\", \"dp_dy_hat\"" << std::endl;
    file << "ZONE T=\"Zone 1\", ZONETYPE=FEQUADRILATERAL, NODES=" << nodes_.size() << ", ELEMENTS=" << n_cells_ << ", DATAPACKING=BLOCK, VARLOCATION=([1,2]=nodal,[3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]=cellcentered), SOLUTIONTIME=" << time << std::endl;

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
        file << std::setw(12) << cells_[i].u_.magnitude() / std::sqrt(cells_[i].p_ / std::pow(cells_[i].gamma_ * cells_[i].p_/std::pow(cells_[i].a_, 2), cells_[i].gamma_)) << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].p_/(cells_[i].gamma_ * cells_[i].p_/std::pow(cells_[i].a_, 2) * R) << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].a_derivative_.x() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].a_derivative_.y() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].ux_derivative_.x() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].ux_derivative_.y() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].uy_derivative_.x() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].uy_derivative_.y() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].p_derivative_.x() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].p_derivative_.y() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].a_derivative_hat_.x() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].a_derivative_hat_.y() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].ux_derivative_hat_.x() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].ux_derivative_hat_.y() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].uy_derivative_hat_.x() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].uy_derivative_hat_.y() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].p_derivative_hat_.x() << " ";
    }
    file << std::endl;

    for (int i = 0; i < n_cells_; ++i) {
        file << std::setw(12) << cells_[i].p_derivative_hat_.y() << " ";
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
        exit(7);
    }

    do {
        std::getline(meshfile, line);  
    }
    while (line.empty());
    
    std::istringstream liness(line);
    liness >> token;
    liness >> value;
    if (token != "NDIME=") {
        std::cerr << "Error: first token should be 'NDIME=', found '" << token << "'. Exiting." << std::endl;
        exit(8);
    }

    if (value != 2) {
        std::cerr << "Error: program only works for 2 dimensions, found '" << value << "'. Exiting." << std::endl;
        exit(9);
    }

    std::vector<Cell_t> farfield;
    std::vector<Cell_t> wall;

    while (!meshfile.eof()) {
        do {
            std::getline(meshfile, line);  
        }
        while (line.empty() && !meshfile.eof());

        std::istringstream liness(line);
        liness >> token;
        std::transform(token.begin(), token.end(), token.begin(),
            [](unsigned char c){ return std::toupper(c); });

        if (token == "NPOIN=") {
            liness >> value;
            nodes_ = std::vector<Node_t>(value);

            for (size_t i = 0; i < nodes_.size(); ++i) {
                std::getline(meshfile, line);
                std::istringstream liness2(line);
                liness2 >> nodes_[i].pos_[0] >> nodes_[i].pos_[1];
            }
        }
        else if (token == "NELEM=") {
            liness >> value;
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
                    exit(10);
                }
            }
        }
        else if (token == "NMARK=") {
            int n_markers;
            liness >> n_markers;

            n_farfield_ = 0;
            n_wall_ = 0;

            for (int i = 0; i < n_markers; ++i) {
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
                std::transform(type.begin(), type.end(), type.begin(),
                    [](unsigned char c){ return std::tolower(c); });

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

                    n_farfield_ += value;
                    farfield.reserve(n_farfield_);

                    for (size_t j = 0; j < value; ++j) {

                        std::getline(meshfile, line);
                        std::istringstream liness6(line);

                        liness6 >> token;
                        if (token != "3") {
                            std::cerr << "Error: expected token '3', found '" << token << "'. Exiting." << std::endl;
                            exit(11);
                        }

                        size_t val0, val1;
                        liness6 >> val0 >> val1;
                        farfield.push_back(Cell_t(2));
                        farfield[farfield.size() - 1].nodes_[0] = val0 - 1;
                        farfield[farfield.size() - 1].nodes_[1] = val1 - 1;
                    }
                }
                else if (type == "wall") {
                    do {
                        std::getline(meshfile, line);
                        if (!line.empty()) {
                            std::istringstream liness(line);
                            liness >> token;
                            liness >> value;
                        }   
                    }
                    while (token != "MARKER_ELEMS=");

                    n_wall_ += value;
                    wall.reserve(n_wall_);

                    for (size_t j = 0; j < value; ++j) {

                        std::getline(meshfile, line);
                        std::istringstream liness6(line);

                        liness6 >> token;
                        if (token != "3") {
                            std::cerr << "Error: expected token '3', found '" << token << "'. Exiting." << std::endl;
                            exit(12);
                        }

                        size_t val0, val1;
                        liness6 >> val0 >> val1;
                        wall.push_back(Cell_t(2));
                        wall[wall.size() - 1].nodes_[0] = val0 - 1;
                        wall[wall.size() - 1].nodes_[1] = val1 - 1;
                    }
                }
                else {
                    std::cerr << "Error: expected marker tag 'farfield' or 'wall', found '" << type << "'. Exiting." << std::endl;
                    exit(6);
                }
            }
        }
        else {
            if (!meshfile.eof()) {
                std::cerr << "Error: expected marker 'NPOIN=', 'NELEM=' or 'NMARK=', found '" << token << "'. Exiting." << std::endl;
                exit(13);
            }
        }
    }

    meshfile.close();

    cells_.insert(std::end(cells_), std::begin(farfield), std::end(farfield));
    cells_.insert(std::end(cells_), std::begin(wall), std::end(wall));
}

void FVM::Entities::Mesh2D_t::reconstruction() {
    #pragma omp parallel for schedule(static)
    for (long long i = 0; i < n_cells_; ++i) {
        FVM::Entities::Cell_t& cell = cells_[i];

        double delta_x2 = 0;
        double delta_y2 = 0;
        double delta_xy = 0;

        double delta_ax = 0;
        double delta_ay = 0;
        double delta_uxx = 0;
        double delta_uxy = 0;
        double delta_uyx = 0;
        double delta_uyy = 0;
        double delta_px = 0;
        double delta_py = 0;
        double delta_gammax = 0; // CHECK how to do those?
        double delta_gammay = 0;
        
        for (auto cell_index: cell.cells_) {
            FVM::Entities::Cell_t& cell_k = cells_[cell_index];

            const double delta_x = cell_k.center_.x() - cell.center_.x();
            const double delta_y = cell_k.center_.y() - cell.center_.y();

            delta_x2 += std::pow(delta_x, 2);
            delta_y2 += std::pow(delta_y, 2);
            delta_xy += delta_x * delta_y;

            delta_ax += (cell_k.a_ - cell.a_) * delta_x;
            delta_ay += (cell_k.a_ - cell.a_) * delta_y;
            delta_uxx += (cell_k.u_.x() - cell.u_.x()) * delta_x;
            delta_uxy += (cell_k.u_.x() - cell.u_.x()) * delta_y;
            delta_uyx += (cell_k.u_.y() - cell.u_.y()) * delta_x;
            delta_uyy += (cell_k.u_.y() - cell.u_.y()) * delta_y;
            delta_px += (cell_k.p_ - cell.p_) * delta_x;
            delta_py += (cell_k.p_ - cell.p_) * delta_y;
            delta_gammax += (cell_k.gamma_ - cell.gamma_) * delta_x;
            delta_gammay += (cell_k.gamma_ - cell.gamma_) * delta_y;
        }

        cell.a_derivative_[1] = delta_ay/(std::pow(delta_xy, 2)/delta_x2 - delta_y2) - delta_ax * delta_xy/(delta_x2 * (std::pow(delta_xy, 2)/delta_x2 - delta_y2));
        cell.a_derivative_[0] = -delta_ax/delta_x2 - cell.a_derivative_[1] * delta_xy/delta_x2;
        cell.ux_derivative_[1] = delta_uxy/(std::pow(delta_xy, 2)/delta_x2 - delta_y2) - delta_uxx * delta_xy/(delta_x2 * (std::pow(delta_xy, 2)/delta_x2 - delta_y2));
        cell.ux_derivative_[0] = -delta_uxx/delta_x2 - cell.ux_derivative_[1] * delta_xy/delta_x2;
        cell.uy_derivative_[1] = delta_uyy/(std::pow(delta_xy, 2)/delta_x2 - delta_y2) - delta_uyx * delta_xy/(delta_x2 * (std::pow(delta_xy, 2)/delta_x2 - delta_y2));
        cell.uy_derivative_[0] = -delta_uyx/delta_x2 - cell.uy_derivative_[1] * delta_xy/delta_x2;
        cell.p_derivative_[1] = delta_py/(std::pow(delta_xy, 2)/delta_x2 - delta_y2) - delta_px * delta_xy/(delta_x2 * (std::pow(delta_xy, 2)/delta_x2 - delta_y2));
        cell.p_derivative_[0] = -delta_px/delta_x2 - cell.p_derivative_[1] * delta_xy/delta_x2;
        cell.gamma_derivative_[1] = delta_gammay/(std::pow(delta_xy, 2)/delta_x2 - delta_y2) - delta_gammax * delta_xy/(delta_x2 * (std::pow(delta_xy, 2)/delta_x2 - delta_y2));
        cell.gamma_derivative_[0] = -delta_gammax/delta_x2 - cell.gamma_derivative_[1] * delta_xy/delta_x2;

        cell.a_derivative_ *= -1;
        cell.ux_derivative_ *= -1;
        cell.uy_derivative_ *= -1;
        cell.p_derivative_ *= -1;
        cell.gamma_derivative_ *= -1;
    }
}

void FVM::Entities::Mesh2D_t::reconstruction_hat() {
    #pragma omp parallel for schedule(static)
    for (long long i = 0; i < n_cells_; ++i) {
        FVM::Entities::Cell_t& cell = cells_[i];

        double delta_x2 = 0;
        double delta_y2 = 0;
        double delta_xy = 0;

        double delta_ax = 0;
        double delta_ay = 0;
        double delta_uxx = 0;
        double delta_uxy = 0;
        double delta_uyx = 0;
        double delta_uyy = 0;
        double delta_px = 0;
        double delta_py = 0;
        double delta_gammax = 0; // CHECK how to do those?
        double delta_gammay = 0;
        
        for (auto cell_index: cell.cells_) {
            FVM::Entities::Cell_t& cell_k = cells_[cell_index];

            const double delta_x = cell_k.center_.x() - cell.center_.x();
            const double delta_y = cell_k.center_.y() - cell.center_.y();

            delta_x2 += std::pow(delta_x, 2);
            delta_y2 += std::pow(delta_y, 2);
            delta_xy += delta_x * delta_y;

            delta_ax += (cell_k.a_hat_ - cell.a_hat_) * delta_x;
            delta_ay += (cell_k.a_hat_ - cell.a_hat_) * delta_y;
            delta_uxx += (cell_k.u_hat_.x() - cell.u_hat_.x()) * delta_x;
            delta_uxy += (cell_k.u_hat_.x() - cell.u_hat_.x()) * delta_y;
            delta_uyx += (cell_k.u_hat_.y() - cell.u_hat_.y()) * delta_x;
            delta_uyy += (cell_k.u_hat_.y() - cell.u_hat_.y()) * delta_y;
            delta_px += (cell_k.p_hat_ - cell.p_hat_) * delta_x;
            delta_py += (cell_k.p_hat_ - cell.p_hat_) * delta_y;
            delta_gammax += (cell_k.gamma_hat_ - cell.gamma_hat_) * delta_x;
            delta_gammay += (cell_k.gamma_hat_ - cell.gamma_hat_) * delta_y;
        }

        cell.a_derivative_hat_[1] = delta_ay/(std::pow(delta_xy, 2)/delta_x2 - delta_y2) - delta_ax * delta_xy/(delta_x2 * (std::pow(delta_xy, 2)/delta_x2 - delta_y2));
        cell.a_derivative_hat_[0] = -delta_ax/delta_x2 - cell.a_derivative_hat_[1] * delta_xy/delta_x2;
        cell.ux_derivative_hat_[1] = delta_uxy/(std::pow(delta_xy, 2)/delta_x2 - delta_y2) - delta_uxx * delta_xy/(delta_x2 * (std::pow(delta_xy, 2)/delta_x2 - delta_y2));
        cell.ux_derivative_hat_[0] = -delta_uxx/delta_x2 - cell.ux_derivative_hat_[1] * delta_xy/delta_x2;
        cell.uy_derivative_hat_[1] = delta_uyy/(std::pow(delta_xy, 2)/delta_x2 - delta_y2) - delta_uyx * delta_xy/(delta_x2 * (std::pow(delta_xy, 2)/delta_x2 - delta_y2));
        cell.uy_derivative_hat_[0] = -delta_uyx/delta_x2 - cell.uy_derivative_hat_[1] * delta_xy/delta_x2;
        cell.p_derivative_hat_[1] = delta_py/(std::pow(delta_xy, 2)/delta_x2 - delta_y2) - delta_px * delta_xy/(delta_x2 * (std::pow(delta_xy, 2)/delta_x2 - delta_y2));
        cell.p_derivative_hat_[0] = -delta_px/delta_x2 - cell.p_derivative_hat_[1] * delta_xy/delta_x2;
        cell.gamma_derivative_hat_[1] = delta_gammay/(std::pow(delta_xy, 2)/delta_x2 - delta_y2) - delta_gammax * delta_xy/(delta_x2 * (std::pow(delta_xy, 2)/delta_x2 - delta_y2));
        cell.gamma_derivative_hat_[0] = -delta_gammax/delta_x2 - cell.gamma_derivative_hat_[1] * delta_xy/delta_x2;

        cell.a_derivative_hat_ *= -1;
        cell.ux_derivative_hat_ *= -1;
        cell.uy_derivative_hat_ *= -1;
        cell.p_derivative_hat_ *= -1;
        cell.gamma_derivative_hat_ *= -1;
    }
}

void FVM::Entities::Mesh2D_t::build_node_to_cell() {
    for (size_t j = 0; j < cells_.size(); ++j) {
        for (auto node_index: cells_[j].nodes_) {
            nodes_[node_index].cells_.push_back(j);
        }
    }
}

void FVM::Entities::Mesh2D_t::build_cell_to_cell() {
    #pragma omp parallel for schedule(static)
    for (long long i = 0; i < cells_.size(); ++i) {
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
        for (size_t j = 0; j < cells_[i].nodes_.size(); ++j) {
            size_t nodes[] = {cells_[i].nodes_[j], (j < cells_[i].nodes_.size() - 1) ? cells_[i].nodes_[j + 1] : cells_[i].nodes_[0]};
            bool found = false;
            for (auto face_index: nodes_[nodes[0]].faces_) {
                if (faces_[face_index].nodes_[0] == nodes[1] || faces_[face_index].nodes_[1] == nodes[1] ) {
                    found = true;
                    faces_[face_index].cells_[1] = i;
                    cells_[i].faces_[j] = face_index;
                    break;
                }
            }

            if (!found) {
                cells_[i].faces_[j] = faces_.size();
                nodes_[nodes[0]].faces_.push_back(faces_.size());
                nodes_[nodes[1]].faces_.push_back(faces_.size());
                faces_.push_back(Face_t(nodes[0], nodes[1], i, -1));
            }
        }
    }
}

void FVM::Entities::Mesh2D_t::compute_cell_geometry() {
    #pragma omp parallel for schedule(static)
    for (long long i = 0; i < cells_.size(); ++i) {
        FVM::Entities::Cell_t& cell = cells_[i];
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
    #pragma omp parallel for schedule(static)
    for (long long i = 0; i < faces_.size(); ++i) {
        FVM::Entities::Face_t& face = faces_[i];
        face.tangent_ = nodes_[face.nodes_[1]].pos_ - nodes_[face.nodes_[0]].pos_; 
        face.length_ = face.tangent_.magnitude();
        face.tangent_ /= face.length_; // CHECK should be normalized or not?
        face.normal_ = Vec2f(face.tangent_.y(), -face.tangent_.x());         

        face.center_ = (nodes_[face.nodes_[0]].pos_ + nodes_[face.nodes_[1]].pos_) * 0.5;
        const Vec2f delta = face.center_ - cells_[face.cells_[0]].center_; // CHECK doesn't work with ghost cells
        const double sign = std::copysign(1.0, face.normal_.dot(delta));
        face.normal_ *= sign;
        face.tangent_ *= sign;
    }
}
