#ifndef FVM_NODE_T_H
#define FVM_NODE_T_H

#include "entities/Vec2f.h"
#include <vector>

namespace FVM { namespace Entities {
    class Node_t { 
        public: 
            Node_t();
            ~Node_t();

            FVM::Entities::Vec2f pos_;
            std::vector<size_t> cells_;
    };
}}
#endif