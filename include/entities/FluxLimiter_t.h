#ifndef FVM_FLUXLIMITER_T_H
#define FVM_FLUXLIMITER_T_H

#include <vector>

namespace FVM { namespace Entities {
    class FluxLimiter_t { 
        public: 
            FluxLimiter_t() {};
            virtual ~FluxLimiter_t() {};

            virtual void calculate_derivatives(const std::vector<double> &x, const std::vector<double> &u, const std::vector<double> &a, const std::vector<double> &p, std::vector<double> &du_dx, std::vector<double> &da_dx, std::vector<double> &dp_dx) = 0;
    };
}}
#endif