#ifndef __SCALAR_MATH_H_
#define __SCALAR_MATH_H_

/**
 * \file common/scalar_math.h
 * \brief Scalar math functions to be used in simulation code.
 *
 * To avoid name collisions with globally-defined functions (as CUDA does), the
 * functions here are suffixed with `r` for ::real. It is strongly recommended
 * to use these functions, because (a) they correctly point to double/float
 * whatever the choice of ::USE_DOUBLE and (b) they will use the correct CPU/GPU
 * version.
 */

#include <cmath>
#include <algorithm>

/// Clamp an integer
inline PHYSICS int clampi(int i, int low, int high)
{
	return i < low ? low : (i > high ? high : i);
}
/// Clamp a real number
inline PHYSICS real clampr(real i, real low, real high)
{
	return i < low ? low : (i > high ? high : i);
}

/// Round to nearest integer, rounding halfway cases to the nearest even number
inline PHYSICS int rintr(real i)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return __double2int_rn(i);
	#else // USE_DOUBLE
		return __float2int_rn(i);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return (int)std::rint(i);
#endif // CUDA_COMPILING
}

/// Round to down to integer (towards negative infinity)
inline PHYSICS int floorr(real i)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return __double2int_rd(i);
	#else // USE_DOUBLE
		return __float2int_rd(i);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return (int)std::floor(i);
#endif // CUDA_COMPILING
}

/// Get the lowest of two real numbers
inline PHYSICS real minr(real a, real b)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return fmin(a, b);
	#else // USE_DOUBLE
		return fminf(a, b);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::min(a, b);
#endif // CUDA_COMPILING
}

/// Get the highest of two real numbers
inline PHYSICS real maxr(real a, real b)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return fmax(a, b);
	#else // USE_DOUBLE
		return fmaxf(a, b);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::max(a, b);
#endif // CUDA_COMPILING
}

/// Take the absolute value of a real number
inline PHYSICS real absr(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return fabs(x);
	#else // USE_DOUBLE
		return fabsf(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::abs(x);
#endif // CUDA_COMPILING
}

/// Take the square root of a real number
inline PHYSICS real sqrtr(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return sqrt(x);
	#else // USE_DOUBLE
		return sqrtf(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::sqrt(x);
#endif // CUDA_COMPILING
}

/// Take the cube root of a real number, even if that number is negative.
inline PHYSICS real cbrtr(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return cbrt(x);
	#else // USE_DOUBLE
		return cbrtf(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::cbrt(x);
#endif // CUDA_COMPILING
}

/// Take the natural logarithm of a real number
inline PHYSICS real logr(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return log(x);
	#else // USE_DOUBLE
		return logf(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::log(x);
#endif // CUDA_COMPILING
}

/// Take the exponent of a real number
inline PHYSICS real expr(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return exp(x);
	#else // USE_DOUBLE
		return expf(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::exp(x);
#endif // CUDA_COMPILING
}

/// Take the exponent of a real number, minus 1: (e^x) - 1
inline PHYSICS real expm1r(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return expm1(x);
	#else // USE_DOUBLE
		return expm1f(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::expm1(x);
#endif // CUDA_COMPILING
}

/// Take the sine of a real number
inline PHYSICS real sinr(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return sin(x);
	#else // USE_DOUBLE
		return sinf(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::sin(x);
#endif // CUDA_COMPILING
}

/// Take the cosine of a real number
inline PHYSICS real cosr(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return cos(x);
	#else // USE_DOUBLE
		return cosf(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::cos(x);
#endif // CUDA_COMPILING
}

/**
 * \brief Calculate the arc tangent of the ratio of the input arguments, in the
 *        correct quadrant as determined by the signs of the parameters.
 */
inline PHYSICS real atan2r(real y, real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return atan2(y, x);
	#else // USE_DOUBLE
		return atan2f(y, x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::atan2(y, x);
#endif // CUDA_COMPILING
}

/**
 * \brief Simultaneously calculate the sine and cosine of a real number.
 *
 * \param x       Number to take the sine and cosine of
 * \param sin_ptr Pointer to the value to store the sine in
 * \param cos_ptr Pointer to the value to store the cosine in
 */
inline PHYSICS void sincosr(real x, real* sin_ptr, real* cos_ptr)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		sincos(x, sin_ptr, cos_ptr);
	#else // USE_DOUBLE
		sincosf(x, sin_ptr, cos_ptr);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	// TODO: try using sin^2(x) + cos^2(x) = 1,
	// but remember to get the sign right.
	*sin_ptr = std::sin(x);
	*cos_ptr = std::cos(x);
#endif // CUDA_COMPILING
}

/**
 * \brief Create value with magnitude, copying sign of another value.
 *
 * Similar to taking sign(y) * abs(x).
 *
 * \param x Number to take magnitude of
 * \param y Number to copy the sign of
 */
inline PHYSICS real copysignr(real x, real y)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return copysign(x, y);
	#else // USE_DOUBLE
		return copysignf(x, y);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	return std::copysign(x, y);
#endif // CUDA_COMPILING
}

/// Clamp a real number between 0 and 1.
inline PHYSICS real saturater(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return clampr(x, 0, 1);
	#else // USE_DOUBLE
		return __saturatef(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	return clampr(x, 0, 1);
#endif // CUDA_COMPILING
}

#endif // __SCALAR_MATH_H_
