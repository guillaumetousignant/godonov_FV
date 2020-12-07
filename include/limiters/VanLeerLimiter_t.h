#ifndef FVM_VANLEERLIMITER_T_H
#define FVM_VANLEERLIMITER_T_H

#include "entities/FluxLimiter_t.h"

namespace FVM { namespace Limiters {
    class VanLeerLimiter_t final : public FVM::Entities::FluxLimiter_t { 
        public: 
            VanLeerLimiter_t();
            virtual ~VanLeerLimiter_t();

            virtual double phi(double a, double b);
    };
}}
#endif