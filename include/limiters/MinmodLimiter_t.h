#ifndef FVM_MINMODLIMITER_T_H
#define FVM_MINMODLIMITER_T_H

#include "entities/FluxLimiter_t.h"

namespace FVM { namespace Limiters {
    class MinmodLimiter_t final : public FVM::Entities::FluxLimiter_t { 
        public: 
            MinmodLimiter_t();
            virtual ~MinmodLimiter_t();

            virtual void phi(double a, double b);
    };
}}
#endif