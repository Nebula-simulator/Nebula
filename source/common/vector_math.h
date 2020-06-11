#ifndef __VECTOR_MATH_H_
#define __VECTOR_MATH_H_

/**
 * \file common/vector_math.h
 * \brief Vector math functions to be used in simulation code.
 */

/// Vector addition
inline PHYSICS vec3 operator+(vec3 a, vec3 b)
{
	return{ a.x + b.x, a.y + b.y, a.z + b.z };
}

/// Vector subtraction
inline PHYSICS vec3 operator-(vec3 a, vec3 b)
{
	return{ a.x - b.x, a.y - b.y, a.z - b.z };
}

/// Vector-scalar multiplication
inline PHYSICS vec3 operator*(vec3 a, real b)
{
	return{ a.x*b, a.y*b, a.z*b };
}

/// Scalar-vector multiplication
inline PHYSICS vec3 operator*(real a, vec3 b)
{
	return{ a*b.x, a*b.y, a*b.z };
}

/// Vector-scalar division
inline PHYSICS vec3 operator/(vec3 a, real b)
{
	return a * (1 / b);
}

/// Negate vector
inline PHYSICS vec3 operator-(vec3 a)
{
	return{ -a.x, -a.y, -a.z };
}

/// Compound vector addition and assignment
inline PHYSICS void operator+=(vec3& a, vec3 b)
{
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
}

/// Compound vector subtraction and assignment
inline PHYSICS void operator-=(vec3& a, vec3 b)
{
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
}

/// Compound vector-scalar multiplication and assignment
inline PHYSICS void operator*=(vec3& a, real b)
{
	a.x *= b;
	a.y *= b;
	a.z *= b;
}

/// Compound vector-scalar division and assignment
inline PHYSICS void operator/=(vec3& a, real b)
{
	a *= 1 / b;
}

/// Vector equality comparison
inline PHYSICS bool operator==(vec3 const & a, vec3 const & b)
{
	return (a.x == b.x
		&& a.y == b.y
		&& a.z == b.z);
}
/// Vector equality comparison
inline PHYSICS bool operator!=(vec3 const & a, vec3 const & b)
{
	return (a.x != b.x
		|| a.y != b.y
		|| a.z != b.z);
}

/// Outer product, or cross product, between vectors
inline PHYSICS vec3 cross_product(vec3 a, vec3 b)
{
	return
	{
		a.y*b.z - a.z*b.y,
		a.z*b.x - a.x*b.z,
		a.x*b.y - a.y*b.x
	};
}

/// Inner product, or dot product, between vectors
inline PHYSICS real dot_product(vec3 a, vec3 b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

/// Magnitude of a vector
inline PHYSICS real magnitude(vec3 a)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return norm3d(a.x, a.y, a.z);
	#else // USE_DOUBLE
		return norm3df(a.x, a.y, a.z);
	#endif //USE_DOUBLE
#else // CUDA_COMPILING
	return sqrtr(a.x*a.x + a.y*a.y + a.z*a.z);
#endif // CUDA_COMPILING
}

/// Magnitude of a vector squared. Faster than taking the magnitude.
inline PHYSICS real magnitude_squared(vec3 a)
{
	return a.x*a.x + a.y*a.y + a.z*a.z;
}

/// Get normalised version of vector, that is, rescaled such that its magnitude is 1.
inline PHYSICS vec3 normalised(vec3 a)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		const real rnorm = rnorm3d(a.x, a.y, a.z);
	#else // USE_DOUBLE
		const real rnorm = rnorm3df(a.x, a.y, a.z);
	#endif //USE_DOUBLE
	return a * rnorm;
#else // CUDA_COMPILING
	return a / sqrtr(a.x*a.x + a.y*a.y + a.z*a.z);
#endif // CUDA_COMPILING
}

/**
 * \brief Normalise the vector in-place
 *
 * `normalise(a)` is equivalent to `a = normalised(a)`.
 */
inline PHYSICS void normalise(vec3& a)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		a *= rnorm3d(a.x, a.y, a.z);
	#else // USE_DOUBLE
		a *= rnorm3df(a.x, a.y, a.z);
	#endif //USE_DOUBLE
#else // CUDA_COMPILING
	a /= sqrtr(a.x*a.x + a.y*a.y + a.z*a.z);
#endif // CUDA_COMPILING
}

/**
 * \brief Find a vector normal to \p dir.
 *
 * This function finds a coordinate system normal to \p dir and returns a vector
 * with angle \p phi in that system. If the input direction is is a unit vector,
 * the returned vector is also a unit vector. Otherwise, the returned vector is
 * not a unit vector, but it is normal to the input direction.
 *
 * \param dir Vector to find a normal vector to.
 * \param phi Rotation angle.
 */
inline PHYSICS vec3 make_normal_vec(vec3 dir, real phi)
{
	real sin_azimuth, cos_azimuth;
	sincosr(atan2r(dir.y, dir.x), &sin_azimuth, &cos_azimuth);

	const vec3 unit_v {
		dir.z*cos_azimuth,
		dir.z*sin_azimuth,
		-sqrtr(dir.x*dir.x + dir.y*dir.y)
	};
	const vec3 unit_u = cross_product(unit_v, dir);

	real sin_phi, cos_phi;
	sincosr(phi, &sin_phi, &cos_phi);
	return unit_u*cos_phi + unit_v*sin_phi;
}

/**
 * \brief Create a unit vector, given Euler angles.
 *
 * \param cos_theta Cosine of the angle w.r.t the z axis.
 * \param phi       Angle in the xy plane.
 */
inline PHYSICS vec3 make_unit_vec(real cos_theta, real phi)
{
	real sin_phi, cos_phi;
	sincosr(phi, &sin_phi, &cos_phi);

	cos_theta = clampr(cos_theta, -1, 1);
	const real sin_theta = sqrtf(1 - cos_theta*cos_theta);
	return vec3 {
		sin_theta * cos_phi,
		sin_theta * sin_phi,
		cos_theta
	};
}

#endif // VECTOR_MATH_H_
