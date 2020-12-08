#ifndef FVM_MINMODLIMITER_T_H
#define FVM_MINMODLIMITER_T_H

#include "entities/FluxLimiter_t.h"

namespace FVM { namespace Limiters {
    class MinmodLimiter_t final : public FVM::Entities::FluxLimiter_t { 
        public: 
            MinmodLimiter_t();
            virtual ~MinmodLimiter_t();

            virtual void calculate_derivatives(const std::vector<double> &x, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &du_dx, std::vector<double> &da_dx, std::vector<double> &dp_dx);
    };
}}
#endif