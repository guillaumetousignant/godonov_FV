#ifndef FVM_ENTITIES_H
#define FVM_ENTITIES_H

namespace FVM { 
    /**
     * @brief Contains the basic entities used through the library.
     * 
     * Most of the entities represent interface classes, to be specialized in their respective namespaces.
     * Some entities, such as Vec3f, represent basic types used throughout the program.
     */
    namespace Entities {}
}

#include "FluxCalculator_t.h"
#include "FluxLimiter_t.h"
#include "Mesh1D_t.h"
#include "Mesh2D_t.h"
#include "Solver_t.h"
#include "Vec2f.h"

#endif