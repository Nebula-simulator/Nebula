#ifndef __GEOMETRY_OCTREE_H_
#define __GEOMETRY_OCTREE_H_

#include "../core/triangle.h"
#include "../core/events.h"
#include "octree/octree_builder.h"

namespace nbl { namespace geometry {

namespace detail
{
	/**
	 * \brief Responsible for managing the memory on the CPU or GPU device.
	 */
	template<bool gpu_flag>
	struct octree_factory;
}

/**
 * \brief Stores geometry as an octree.
 *
 * This class is responsible for the collision detection system. It holds the
 * simulation domain (a finite axis-aligned box) and the triangles.
 *
 * This is an octree implementation, which is much faster for large numbers of
 * triangles than a simple list.
 *
 * It is allocated by the static {@link create} and {@link destroy} functions,
 * there is no constructor or destructor to be used.
 */
template<bool gpu_flag>
class octree
{
public:
	using triangle_index_t = uint32_t; ///< Type for indexing the triangles

	/**
	 * \brief Send the octree to the correct device
	 *
	 * \param linoct Linearized octree to send to the device.
	 */
	static CPU octree<gpu_flag> create(nbl::geometry::octree_builder::linearized_octree const & linoct);

	/**
	 * \brief Build the octree from a list of triangles, and send it to the correct device
	 *
	 * \param triangles List of triangles to be used in the simulation.
	 */
	static CPU octree<gpu_flag> create(std::vector<triangle> const & triangles);

	/**
	 * \brief Destroy the octree, deallocating the data.
	 */
	static CPU void destroy(octree & geometry);

	/**
	 * \brief Check whether a certain position is part of the simulation domain.
	 *
	 * \param position Position to check.
	 */
	inline PHYSICS bool in_domain(vec3 position) const;

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
	inline static PHYSICS vec3 AABB_intersect(vec3 pos, vec3 dir, vec3 center, vec3 halfsize);

	inline static PHYSICS int clz(uint64_t x);

	int* _octree_data    = nullptr;
	triangle* _triangles = nullptr;
	vec3 _AABB_center    = { 0, 0, 0 };
	vec3 _AABB_halfsize  = { 0, 0, 0 };
	real _max_extent     = 0;

	friend struct detail::octree_factory<gpu_flag>;
};

}} // namespace nbl::geometry

#include "octree.inl"

#endif // __GEOMETRY_OCTREE_H_
