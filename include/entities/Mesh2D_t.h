#ifndef FVM_MESH2D_T_H
#define FVM_MESH2D_T_H

#include "entities/Vec2f.h"
#include "entities/Cell_t.h"
#include "entities/Face_t.h"
#include "entities/Node_t.h"
#include <vector>
#include <string>
#include <filesystem>


struct state { 
    double rho;
    FVM::Entities::Vec2f u;
    double p;
    double gamma; 
};

namespace FVM { namespace Entities {
    class Mesh2D_t { 
        public: 
            Mesh2D_t(std::filesystem::path filename);

            size_t n_cells_;
            size_t n_farfield_;
            size_t n_wall_;
            std::vector<Node_t> nodes_; // x, y
            std::vector<Cell_t> cells_;
            std::vector<Face_t> faces_;

            void initial_conditions(FVM::Entities::Vec2f center, const state& state_NE, const state& state_NW, const state& state_SW, const state& state_SE);
            void boundary_conditions();
            void boundary_conditions_hat();
            void write_tecplot(std::filesystem::path filename, int problem_number, double time);
            void reconstruction();
            void reconstruction_hat();

        private:
            void read_su2(std::filesystem::path filename);
            void build_node_to_cell();
            void build_cell_to_cell();
            void build_faces();
            void compute_cell_geometry();
            void compute_face_geometry();
    };
}}
#endif