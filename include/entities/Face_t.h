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
            FVM::Entities::Vec2f center_;

            // State
            double F_1_;
            double F_2_;
            double F_3_;
            double F_4_;

    };
}}
#endif