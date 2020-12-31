#ifndef FVM_VENKATAKRISHNANLIMITER_T_H
#define FVM_VENKATAKRISHNANLIMITER_T_H

#include "entities/Mesh2D_t.h"

namespace FVM { namespace Limiters {
    class VenkatakrishnanLimiter_t { 
        public: 
            VenkatakrishnanLimiter_t();
            ~VenkatakrishnanLimiter_t();

            void calculate_derivatives(FVM::Entities::Mesh2D_t &mesh);
    };
}}
#endif