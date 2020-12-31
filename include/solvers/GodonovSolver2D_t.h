#ifndef FVM_GODONOVSOLVER2D_T_H
#define FVM_GODONOVSOLVER2D_T_H

#include "entities/Mesh2D_t.h"
#include <vector>
#include <string>

namespace FVM { namespace Solvers {
    template<typename FluxCalculator, typename FluxLimiter>
    class GodonovSolver2D_t { 
        public: 
            GodonovSolver2D_t();

            FluxCalculator flux_calculator_;
            FluxLimiter flux_limiter_;

            void solve(double end_time, double cfl, FVM::Entities::Mesh2D_t& mesh);

        private:
            double calculate_delta_t(double cfl, FVM::Entities::Mesh2D_t& mesh);
            void timestep(double delta_t, FVM::Entities::Mesh2D_t& mesh);
    };
}}
#endif