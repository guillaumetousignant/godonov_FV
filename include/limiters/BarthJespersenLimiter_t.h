#ifndef FVM_BARTHJESPERSENLIMITER_T_H
#define FVM_BARTHJESPERSENLIMITER_T_H

#include "entities/Mesh2D_t.h"

namespace FVM { namespace Limiters {
    class BarthJespersenLimiter_t { 
        public: 
            BarthJespersenLimiter_t();
            ~BarthJespersenLimiter_t();

            void calculate_derivatives(FVM::Entities::Mesh2D_t &mesh);
    };
}}
#endif