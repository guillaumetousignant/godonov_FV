#ifndef FVM_ELEMENT_T_H
#define FVM_ELEMENT_T_H

#include <vector>

namespace FVM { namespace Entities {
    class Element_t { 
        public: 
            Element_t(int n_sides);
            ~Element_t();

            std::vector<int> nodes_;
    };
}}
#endif