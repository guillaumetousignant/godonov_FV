#ifndef FVM_FACE_T_H
#define FVM_FACE_T_H

#include <vector>
#include "entities/Vec2f.h"

namespace FVM { namespace Entities {
    class Face_t { 
        public: 
            Face_t();
            Face_t(size_t node_0, size_t node_1, size_t cell_0, size_t cell_1);
            ~Face_t();

            // Connectivity
            size_t nodes_[2];
            size_t cells_[2];

            // Geometry
            FVM::Entities::Vec2f normal_; // Not normalized, norm is the area
            FVM::Entities::Vec2f tangent_; // Not normalized, norm is the area

    };
}}
#endif