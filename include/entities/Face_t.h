#ifndef FVM_FACE_T_H
#define FVM_FACE_T_H

#include <vector>

namespace FVM { namespace Entities {
    class Face_t { 
        public: 
            Face_t();
            Face_t(size_t node_0, size_t node_1, size_t cell_0, size_t cell_1);
            ~Face_t();

            size_t nodes_[2];
            size_t cells_[2];
    };
}}
#endif