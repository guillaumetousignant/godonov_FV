#ifndef FVM_MESH2D_T_H
#define FVM_MESH2D_T_H

#include "entities/Vec2f.h"
#include "entities/Element_t.h"
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
            std::vector<Vec2f> points_;
            std::vector<Element_t> cells_;

            std::vector<double> a_;
            std::vector<double> u_;
            std::vector<double> p_;
            std::vector<double> x_;
            std::vector<double> gamma_;
            double delta_x_;
            int n_faces_;
            std::vector<double> F_1_; // ρu
            std::vector<double> F_2_; // ρu² + p
            std::vector<double> F_3_; // u(γp/(γ-1) + ρu²/2)

            void initial_conditions(double a_L, double a_R, double u_L, double u_R, double p_L, double p_R, double x_L, double x_R, double gamma_L, double gamma_R, double discontinuity);

        private:
            void readSU2(std::filesystem::path filename);
    };
}}
#endif