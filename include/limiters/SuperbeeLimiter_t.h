#ifndef FVM_SUPERBEELIMITER_T_H
#define FVM_SUPERBEELIMITER_T_H

#include "entities/FluxLimiter_t.h"

namespace FVM { namespace Limiters {
    class SuperbeeLimiter_t final : public FVM::Entities::FluxLimiter_t { 
        public: 
            SuperbeeLimiter_t();
            virtual ~SuperbeeLimiter_t();

            virtual void phi(double a, double b);
    };
}}
#endif