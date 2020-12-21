#ifndef FVM_GODONOVSOLVER2D_T_H
#define FVM_GODONOVSOLVER2D_T_H

#include "entities/Mesh2D_t.h"
#include <vector>
#include <string>

namespace FVM { namespace Solvers {
    template<typename FluxCalculator>
    class GodonovSolver2D_t { 
        public: 
            GodonovSolver2D_t();
            ~GodonovSolver2D_t();

            FluxCalculator flux_calculator_;

            void solve(double time, double cfl, FVM::Entities::Mesh2D_t& mesh);

        private:
            double calculate_delta_t();
            void timestep(double delta_t, double delta_x, const std::vector<double> &gamma, std::vector<double> &u, std::vector<double> &a, std::vector<double> &p, const std::vector<double> F_1, const std::vector<double> F_2, const std::vector<double> F_3);
    };
}}
#endif