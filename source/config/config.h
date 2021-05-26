#ifndef __CONFIG_H_
#define __CONFIG_H_

/**
 * \file config/config.h
 */

/**
 * \brief Whether or not to use double precision in simulations.
 *
 * Set to 0 for single-precision and to 1 for double-precision calculations.
 */
#define USE_DOUBLE 0

// Include relevant headers
#include <cstdint>
#include <cstddef>
#include <stdexcept>

#include "device.h"
#include "data_types.h"

#include "../common/scalar_math.h"
#include "../common/vector_math.h"

#if CUDA_COMPILER_AVAILABLE
	#include "../common/cuda/cuda_new.h"
	#include "../common/cuda/cuda_mem_scope.h"
#endif // CUDA_COMPILER_AVAILABLE

/// A frequently-used constant of mathematics
constexpr real pi = 3.1415926535897932_r;

/// Version number
constexpr unsigned int VERSION_MAJOR = 1;
constexpr unsigned int VERSION_MINOR = 0;
constexpr unsigned int VERSION_PATCH = 2;

#endif // __CONFIG_H_
