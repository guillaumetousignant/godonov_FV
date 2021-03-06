#ifndef FVM_SOLVERS_H
#define FVM_SOLVERS_H

namespace FVM { 
    /**
     * @brief Contains the different types of solver objects that can be used through the library.
     * 
     * Solvers are the core of the program, they compute the solution to the given flow problems.
     */
    namespace Solvers {}
}

#include "ExactSolver_t.h"
#include "GodonovSolver_t.h"
#include "GodonovSolverHigherOrder_t.h"
#include "GodonovSolver2D_t.h"

#endif