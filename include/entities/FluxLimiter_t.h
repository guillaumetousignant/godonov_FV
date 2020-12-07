#ifndef FLUXLIMITER_T_H
#define FLUXLIMITER_T_H

#include <vector>

class FluxLimiter_t { 
    public: 
        FluxLimiter_t() {};
        virtual ~FluxLimiter_t() {};

        virtual void phi(double a, double b) = 0;
};

#endif