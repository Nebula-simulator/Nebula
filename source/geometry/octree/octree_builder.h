#ifndef __OCTREE_BUILDER_H_
#define __OCTREE_BUILDER_H_

#include <vector>
#include "../../config/config.h"
#include "../../core/triangle.h"

namespace nbl { namespace geometry { namespace octree_builder {

class octree_node;
class octree_child;
class octree_root;

/**
 * \brief Octree node base class.
 *
 * There are two specializations: a root node, which has no parents, and a
 * general child node, which does.
 *
 * The present class exists mostly to provide a common interface.
 */
class octree_node
{
public:
	/**
	 * \brief Check if the current node is a leaf node.
	 */
	bool is_leaf() const;

	/**
	 * \brief Get this node's depth
	 *
	 * The root node has depth 0, its children have depth 1, and so on.
	 */
	virtual int depth() const = 0;

	/**
	 * \brief Compute the center and half-size.
	 *
	 * This is done by traversing the tree back to the root.
	 *
	 * \return First element is the center coordinate, second is the half-size.
	 */
	virtual std::pair<vec3, vec3> get_center_halfsize() const = 0;

	/**
	 * \brief Compute the current node's Morton index with respec to the root.
	 *
	 * This is done by traversing the tree back to the root.
	 */
	virtual uint64_t morton() const = 0;

	/**
	 * \brief Get a pointer to the child node.
	 *
	 * If {@link is_leaf()} is `true`, this function will return a valid pointer.
	 * Otherwise, it returns `nullptr`.
	 *
	 * \param octant Child's octant number between 0 and 7.
	 */
	octree_child const * get_child(int octant) const;

	/**
	 * \brief Get a reference to the root node.
	 *
	 * This is done by traversing the tree back to the root.
	 */
	virtual octree_root const & get_root() const = 0;

	virtual ~octree_node();

protected:
	/**
	 * \brief Pointer to the first child node.
	 *
	 * If not null, then there are 8 children at successive memory locations.
	 */
	octree_child* _children = nullptr;

	/**
	 * \brief Test for triangle overlap with this node
	 */
	bool overlap(triangle const & tri) const;

	/**
	 * \brief Triangle-box overlap test
	 */
	static bool overlap(triangle const & tri, vec3 cell_center, vec3 cell_halfsize);

	/**
	 * \brief Get center and half-size extent of a child octant.
	 */
	std::pair<vec3, vec3> get_octant_center_halfsize(int octant) const;
	/**
	 * \brief Get center and half-size extent of a child octant, if the present
	 * node's center and half-size are already known.
	 */
	static std::pair<vec3, vec3> get_octant_center_halfsize(
		int octant, vec3 my_center, vec3 my_halfsize);

	friend class octree_child;
	friend class octree_root;

	octree_node() = default;
	octree_node(octree_node const &) = delete;
	octree_node(octree_node&&) = delete;
	octree_node& operator=(octree_node&) = delete;
	octree_node& operator=(octree_node&&) = delete;
};


class octree_child : public octree_node
{
public:
	virtual int depth() const;

	virtual std::pair<vec3, vec3> get_center_halfsize() const;
	virtual uint64_t morton() const;

	virtual octree_root const & get_root() const;

	/**
	 * \brief Get triangle indices in this node.
	 *
	 * They are indices to get_root()->triangles().
	 * Only works if is_leaf().
	 */
	std::vector<size_t> const & triangles() const;

protected:
	/**
	 * \brief Insert a triangle into this node.
	 *
	 * This function should only be called by the parent node. The parent node
	 * knows this cells' center and halfsize, so this function takes them as
	 * parameters.
	 */
	bool insert(size_t tri_index, triangle const & tri_data, vec3 my_center, vec3 my_halfsize);

private:
	octree_node* _parent = nullptr;

	// If this is a leaf, children is nullptr and _triangles is filled.
	// If not a leaf, children is filled and _triangles is empty.
	std::vector<size_t> _triangles;

	// Get the node this octant represents for its parent.
	int my_octant() const;

	// Split the current leaf node into eight children, and redistribute the triangles.
	void split(vec3 my_center, vec3 my_halfsize);

	// Clear the triangles vector and de-allocate the heap memory as much as the
	// implementation of std::vector allows
	void clear_triangles();

	static constexpr int _max_depth = 21;
	static constexpr int _split_count = 16;

	friend class octree_root;
	octree_child() = default;
	octree_child(octree_child const &) = delete;
	octree_child(octree_child&&) = delete;
	octree_child& operator=(octree_child&) = delete;
	octree_child& operator=(octree_child&&) = delete;
};

class octree_root : public octree_node
{
public:
	/**
	 * \brief Constructor, providing the minimum and maximum coordinates of the
	 *        root box.
	 */
	octree_root(vec3 min, vec3 max);

	/**
	 * \brief Constructor, providing a list of triangles.
	 */
	octree_root(std::vector<triangle> const & triangles);

	virtual int depth() const;

	virtual std::pair<vec3, vec3> get_center_halfsize() const;
	virtual uint64_t morton() const;

	virtual octree_root const & get_root() const { return *this; }

	/**
	 * \brief Insert a new triangle into the octree.
	 *
	 * \return Whether the triangle was successfully inserted.
	 */
	bool insert(triangle const & tri);

	/**
	 * \brief Get access to all triangles stored in the octree.
	 */
	std::vector<triangle> const & triangles() const;

private:
	vec3 center   = { 0, 0, 0 };
	vec3 halfsize = { 0, 0, 0 };
	std::vector<triangle> _triangles;

	friend class linearized_octree;
};

class linearized_octree
{
public:
	/**
	 * \brief Construct to invalid state
	 */
	explicit linearized_octree();

	/**
	 * \brief Construct from a hierarchical octree structure.
	 */
	linearized_octree(octree_root const & root);

	/**
	 * \brief Octree data.
	 *
	 * Contains indices, either into itself (if not a leaf node) or into the
	 * triangle_data vector (if a leaf node).
	 *  * value=0 : child does not exist
	 *  * value>0 : non-leaf child with node indices
	 *  * value<0 : leaf child with triangle indices (index -1 means no triangle)
	 */
	std::vector<int> octree_data;

	/**
	 * \brief Triangle data, approximately sorted by Morton ordering.
	 */
	std::vector<triangle> triangle_data;

	/// Coordinate of the domain center
	vec3 center   = { 0, 0, 0 };

	/// The domain's half-extent
	vec3 halfsize = { 0, 0, 0 };

	/// Get the (axis-aligned) simulation domain.
	vec3 AABB_min() const;
	/// Get the (axis-aligned) simulation domain.
	vec3 AABB_max() const;
};

}}} // namespace nbl::geometry::octree_builder

#endif // __OCTREE_BUILDER_H_
