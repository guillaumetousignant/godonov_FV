#ifndef FVM_CELL_T_H
#define FVM_CELL_T_H

#include <vector>
#include "entities/Vec2f.h"

namespace FVM { namespace Entities {
    class Cell_t { 
        public: 
            Cell_t(); // For use in vectors
            Cell_t(int n_sides);

            // Connectivity
            std::vector<size_t> nodes_;
            std::vector<size_t> cells_;
            std::vector<size_t> faces_;

            // Geometry
            FVM::Entities::Vec2f center_;
            double area_;

            // State
            double a_;
            FVM::Entities::Vec2f u_;
            double p_;
            double gamma_;
            FVM::Entities::Vec2f a_derivative_;
            FVM::Entities::Vec2f ux_derivative_;
            FVM::Entities::Vec2f uy_derivative_;
            FVM::Entities::Vec2f p_derivative_;
            FVM::Entities::Vec2f gamma_derivative_;
            double phi_a_;
            FVM::Entities::Vec2f phi_u_;
            double phi_p_;
            double phi_gamma_;
    };
}}
#endif