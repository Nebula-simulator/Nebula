#ifndef __GEOMETRY_TRILIST_H_
#define __GEOMETRY_TRILIST_H_

#include "../core/triangle.h"
#include "../core/events.h"

namespace nbl { namespace geometry {

namespace detail
{
	/**
	 * \brief Responsible for managing the memory on the CPU or GPU device.
	 */
	template<bool gpu_flag>
	struct trilist_factory;
}

/**
 * \brief Stores geometry as a list of triangles.
 *
 * This class is responsible for the collision detection system. It holds the
 * simulation domain (a finite axis-aligned box) and the triangles.
 *
 * This is the simplest implementation, simply holding a list of all triangles.
 * When checking for collisions, it simply tries all of them in turn.
 *
 * It is allocated by the static {@link create} and {@link destroy} functions,
 * there is no constructor or destructor to be used.
 */
template<bool gpu_flag>
class trilist
{
public:
	using triangle_index_t = uint32_t; ///< Type for indexing the triangles

	/**
	 * \brief Allocate memory for the triangles on the correct device.
	 *
	 * \param triangles List of triangles to be used in the simulation.
	 */
	static CPU trilist create(std::vector<triangle> const & triangles);
	/**
	 * \brief Destroy the triangle list, deallocating the data.
	 */
	static CPU void destroy(trilist & geometry);

	/**
	 * \brief Check whether a certain position is part of the simulation domain.
	 *
	 * \param position Position to check.
	 */
	inline PHYSICS bool in_domain(vec3 position);

	/**
	 * \brief Try to move a particle, checking for collisions with triangles.
	 *
	 * Try to propagate from \p start to \p start + \p distance * \p direction.
	 * If there is no collision with a triangle, the returned intersect_event
	 * contains a `nullptr` triangle pointer.
	 *
	 * \param start           The particle's starting position
	 * \param direction       Direction the particle is going in
	 * \param distance        Distance the particle travels, in units of direction
	 * \param ignore_triangle Triangle to ignore. Used to prevent the particle
	 *                        from intersecting the same triangle twice in case
	 *                        of numerical error.
	 * \param ignore_material Destination material to ignore (most likely the
	 *                        particle's current material).
	 *
	 * TODO ignore_material datatype
	 */
	inline PHYSICS intersect_event propagate(vec3 start, vec3 direction, real distance,
		triangle const * ignore_triangle, int ignore_material) const;

	/**
	 * \brief Get the maximum distance that can be travelled inside the
	 *        simulation domain.
	 */
	inline PHYSICS real get_max_extent() const;

	/**
	 * \brief Get the (axis-aligned) simulation domain.
	 */
	inline PHYSICS vec3 AABB_min() const;
	/**
	 * \brief Get the (axis-aligned) simulation domain.
	 */
	inline PHYSICS vec3 AABB_max() const;

private:
	CPU void set_AABB(vec3 min, vec3 max);

	triangle* _triangles = nullptr;
	triangle_index_t _N  = 0;
	vec3 _AABB_min       = { 0, 0, 0 };
	vec3 _AABB_max       = { 0, 0, 0 };
	real _max_extent     = 0;

	friend struct detail::trilist_factory<gpu_flag>;
};

}} // namespace nbl::geometry

#include "trilist.inl"

#endif // __GEOMETRY_TRILIST_H_
