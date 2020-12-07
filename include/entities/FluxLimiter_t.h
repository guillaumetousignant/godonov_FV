#ifndef FVM_FLUXLIMITER_T_H
#define FVM_FLUXLIMITER_T_H

#include <vector>

namespace FVM { namespace Entities {
    class FluxLimiter_t { 
        public: 
            FluxLimiter_t() {};
            virtual ~FluxLimiter_t() {};

            virtual double phi(double a, double b) = 0;
    };
}}
#endif