#ifndef FVM_H
#define FVM_H

/**
 * @brief Contains all the objects and functions of the library.
 * 
 * All the different namespaces are pulled under this one.
 */
namespace FVM {}

#include "entities/entities.h"
#include "fluxes/fluxes.h"
#include "functions/functions.h"
#include "limiters/limiters.h"
#include "solvers/solvers.h"

#endif