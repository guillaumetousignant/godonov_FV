#ifndef FVM_FLUXES_H
#define FVM_FLUXES_H

namespace FVM { 
    /**
     * @brief Contains the different types of flux calculator objects that can be used through the library.
     * 
     * Flux calculators compute fluxes between cells.
     */
    namespace Fluxes {}
}

#include "ExactRiemannFlux_t.h"
#include "HLLEFlux_t.h"
#include "RoeEntropyFlux_t.h"
#include "RoeFlux_t.h"

#endif