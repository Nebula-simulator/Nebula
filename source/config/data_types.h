#ifndef __DATA_TYPES_H_
#define __DATA_TYPES_H_

/**
 * \file config/data_types.h
 * \brief Define data types to be used in simulation code.
 */

#include <limits>

// Index type
using std::size_t;

/**
 * \typedef real
 * \brief   "Real number" type used in simulations.
 *
 * Define the ::USE_DOUBLE macro as 1 to use double precision,
 * define it as 0 to use single precision.
 */

/**
 * \struct vec3
 * \brief  3-dimensional vector, datatype is ::real.
 */

/**
 * \struct vec4
 * \brief  4-dimensional vector, datatype is ::real.
 */
// scalar and vector types
#if CUDA_HEADERS_AVAILABLE
	#if USE_DOUBLE
		using real = double;
		using vec3 = double3;
		using vec4 = double4;
	#else // USE_DOUBLE
		using real = float;
		using vec3 = float3;
		using vec4 = float4;
	#endif // USE_DOUBLE
#else // CUDA_HEADERS_AVAILABLE
	#if USE_DOUBLE
		using real = double;
		struct vec3 { double x, y, z; };
		struct vec4 { double x, y, z, w; };
	#else // USE_DOUBLE
		using real = float;
		struct vec3 { float x, y, z; };
		struct vec4 { float x, y, z, w; };
	#endif // USE_DOUBLE
#endif // CUDA_HEADERS_AVAILABLE

/// Small value limiting accuracy of numerical computations
constexpr real EPSILON = 10 * std::numeric_limits<real>::epsilon();

/**
 * \brief User-defined literal for using real numbers
 *
 * Example: `1.32_r` is `(real)1.32`.
 */
constexpr PHYSICS real operator"" _r(long double r)
{
	return static_cast<real>(r);
}

#endif // __DATA_TYPES_H_
