#ifndef FVM_FACE_T_H
#define FVM_FACE_T_H

#include <vector>
#include "entities/Vec2f.h"

namespace FVM { namespace Entities {
    class Face_t { 
        public: 
            Face_t();
            Face_t(size_t node_0, size_t node_1, size_t cell_0, size_t cell_1);

            // Connectivity
            size_t nodes_[2];
            size_t cells_[2];

            // Geometry
            FVM::Entities::Vec2f normal_;
            FVM::Entities::Vec2f tangent_;
            double length_;

            // State
            FVM::Entities::Vec2f F_1_;
            FVM::Entities::Vec2f F_2_;
            FVM::Entities::Vec2f F_3_;
            FVM::Entities::Vec2f F_4_;

    };
}}
#endif