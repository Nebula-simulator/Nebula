#ifndef __TRIANGLE_H_
#define __TRIANGLE_H_

/**
 * \brief Triangle class, to be used in simulations.
 */

class triangle
{
public:
	/**
	 * \brief Constructor. takes vertices for the edges.
	 *
	 * TODO: document definition for "in" and "out" materials
	 */
	inline CPU triangle(vec3 r0, vec3 r1, vec3 r2, int material_in, int material_out);

	/**
	 * \brief Find next intersection of a ray with the present triangle.
	 *
	 * The intersection is at `ray_start + t*ray_direction`, where `t` is the
	 * return value of this function.
	 *
	 * Returns -1 if there is no intersection, or if `t` is negative.
	 */
	inline PHYSICS real intersect_ray(vec3 ray_start, vec3 ray_direction) const;

	/**
	 * \brief Get `r0` vertex supplied to constructor.
	 */
	inline PHYSICS vec3 r0() const;
	/**
	 * \brief Get `r1` vertex supplied to constructor.
	 *
	 * This is computed from r0 and the edge between r0 and r1.
	 */
	inline PHYSICS vec3 r1() const;
	/**
	 * \brief Get `r2` vertex supplied to constructor.
	 *
	 * This is computed from r0 and the edge between r0 and r2.
	 */
	inline PHYSICS vec3 r2() const;

	/// Compute the minimal coordinates contained in this triangle.
	inline PHYSICS vec3 AABB_min() const;
	/// Compute the maximal coordinates contained in this triangle.
	inline PHYSICS vec3 AABB_max() const;

	/**
	 * \brief Compute normal vector (unnormalised).
	 *
	 * This is given by `(r1-r0) x (r2-r0)`.
	 */
	inline PHYSICS vec3 get_normal() const;

	/// Equality comparison
	inline PHYSICS bool operator==(triangle const & rhs) const;

	/**
	 * \brief Test for perfect overlap with another triangle.
	 *
	 * This is different from testing equality, for two reasons:
	 *  1. It checks all permutations of the three vertices
	 *  2. It does not check if the materials are the same.
	 */
	inline PHYSICS bool perfect_overlap(triangle const & other) const;

	// TODO: typename material_manager::material_index_t
	int material_in;  ///< Material on side the normal points away from
	int material_out; ///< Material on side the normal points towards

private:
	vec3 _r0;
	vec3 _e1;
	vec3 _e2;
};

#include "triangle.inl"

#endif // __TRIANGLE_H_
