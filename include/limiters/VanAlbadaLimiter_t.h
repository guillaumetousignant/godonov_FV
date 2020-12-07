#ifndef FVM_VANALBADALIMITER_T_H
#define FVM_VANALBADALIMITER_T_H

#include "entities/FluxLimiter_t.h"

namespace FVM { namespace Limiters {
    class VanAlbadaLimiter_t final : public FVM::Entities::FluxLimiter_t { 
        public: 
            VanAlbadaLimiter_t();
            virtual ~VanAlbadaLimiter_t();

            virtual void phi(double a, double b);
    };
}}
#endif