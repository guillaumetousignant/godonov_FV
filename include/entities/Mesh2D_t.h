#ifndef FVM_MESH2D_T_H
#define FVM_MESH2D_T_H

#include "entities/Vec2f.h"
#include "entities/Cell_t.h"
#include "entities/Face_t.h"
#include <vector>
#include <string>
#include <filesystem>

namespace FVM { namespace Entities {
    class Mesh2D_t { 
        public: 
            Mesh2D_t(std::filesystem::path filename);
            ~Mesh2D_t();

            size_t n_cells_;
            size_t n_boundary_;
            std::vector<Vec2f> nodes_; // x, y
            std::vector<Cell_t> cells_;
            std::vector<std::vector<size_t>> node_to_cell_;
            std::vector<Face_t> faces_;

            void initial_conditions(double a_L, double a_R, double u_L, double u_R, double p_L, double p_R, double x_L, double x_R, double gamma_L, double gamma_R, double discontinuity);

        private:
            void readSU2(std::filesystem::path filename);
            void build_node_to_cell();
            void build_cell_to_cell();
            void build_faces();
    };
}}
#endif