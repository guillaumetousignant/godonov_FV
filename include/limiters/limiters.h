#ifndef LIMITERS_H
#define LIMITERS_H

namespace FVM { 
    /**
     * @brief Contains the different types of limiter objects that can be used through the library.
     * 
     * Limiters prevent the apparition of new maxima and minima in the solution.
     */
    namespace Limiters {}
}

#include "MinmodLimiter_t.h"
#include "SuperbeeLimiter_t.h"
#include "VanAlbadaLimiter_t.h"
#include "VanLeerLimiter_t.h"

#endif