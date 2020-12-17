#ifndef FVM_CELL_T_H
#define FVM_CELL_T_H

#include <vector>

namespace FVM { namespace Entities {
    class Cell_t { 
        public: 
            Cell_t(); // For use in vectors
            Cell_t(int n_sides);
            ~Cell_t();

            std::vector<size_t> nodes_;
            std::vector<size_t> cells_;
    };
}}
#endif