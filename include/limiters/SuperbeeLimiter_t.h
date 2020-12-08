#ifndef FVM_SUPERBEELIMITER_T_H
#define FVM_SUPERBEELIMITER_T_H

#include "entities/FluxLimiter_t.h"

namespace FVM { namespace Limiters {
    class SuperbeeLimiter_t final : public FVM::Entities::FluxLimiter_t { 
        public: 
            SuperbeeLimiter_t();
            virtual ~SuperbeeLimiter_t();

            virtual void calculate_derivatives(const std::vector<double> &x, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &du_dx, std::vector<double> &da_dx, std::vector<double> &dp_dx);
    };
}}
#endif