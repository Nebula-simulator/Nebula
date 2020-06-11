#include "octree_builder.h"
#include <assert.h>
#include <queue>
#include <map>
#include <iostream>
#include "tribox.hh"

namespace nbl { namespace geometry { namespace octree_builder {

bool octree_node::is_leaf() const
{
	return (_children == nullptr);
}

octree_child const * octree_node::get_child(int octant) const
{
	assert(octant >= 0 && octant < 8);

	if (_children == nullptr)
		return nullptr;
	return _children + octant;
}

octree_node::~octree_node()
{
	if (_children)
		delete[] _children;
}

bool octree_node::overlap(triangle const & tri, vec3 cell_center, vec3 cell_halfsize)
{
	double center[3] = { cell_center.x, cell_center.y, cell_center.z };
	double halfsize[3] = { 1.01*cell_halfsize.x, 1.01*cell_halfsize.y, 1.01*cell_halfsize.z };
	double triverts[3][3];
	triverts[0][0] = tri.r0().x; triverts[0][1] = tri.r0().y; triverts[0][2] = tri.r0().z;
	triverts[1][0] = tri.r1().x; triverts[1][1] = tri.r1().y; triverts[1][2] = tri.r1().z;
	triverts[2][0] = tri.r2().x; triverts[2][1] = tri.r2().y; triverts[2][2] = tri.r2().z;
	return (triBoxOverlap(center, halfsize, triverts) > 0);
}

bool octree_node::overlap(triangle const & tri) const
{
	const auto my_ch = get_center_halfsize();
	return overlap(tri, my_ch.first, my_ch.second);
}

std::pair<vec3, vec3> octree_node::get_octant_center_halfsize(int octant, vec3 my_center, vec3 my_halfsize)
{
	assert(octant >= 0 && octant < 8);

	const vec3 octant_halfsize = .5 * my_halfsize;
	const vec3 octant_center
	{
		my_center.x + (((octant & 1) << 1) - 1)*octant_halfsize.x,
		my_center.y + (((octant & 2))      - 1)*octant_halfsize.y,
		my_center.z + (((octant & 4) >> 1) - 1)*octant_halfsize.z
	};

	return std::pair<vec3, vec3>
	{
		octant_center,
		octant_halfsize
	};
}

std::pair<vec3, vec3> octree_node::get_octant_center_halfsize(int octant) const
{
	const auto my_ch = get_center_halfsize();
	return get_octant_center_halfsize(octant, my_ch.first, my_ch.second);
}

int octree_child::depth() const
{
	return 1 + _parent->depth();
}

std::pair<vec3, vec3> octree_child::get_center_halfsize() const
{
	assert(_parent != nullptr);
	return _parent->get_octant_center_halfsize(my_octant());
}

uint64_t octree_child::morton() const
{
	return (_parent->morton() << 3) | my_octant();
}

octree_root const & octree_child::get_root() const
{
	return _parent->get_root();
}

bool octree_child::insert(size_t tri_index,
	triangle const & tri_data,
	vec3 my_center, vec3 my_halfsize)
{
	if (!overlap(tri_data, my_center, my_halfsize))
		return false;

	if (is_leaf())
	{
		std::vector<triangle> const & tris = get_root().triangles();
		for (auto tri_idx : _triangles)
		{
			if (tri_data.perfect_overlap(tris[tri_idx]))
			{
				std::clog << "WARNING: ignoring duplicate triangle!" << std::endl;
				return false;
			}
		}

		_triangles.push_back(tri_index);
		if (_triangles.size() > _split_count && depth() < _max_depth)
			split(my_center, my_halfsize);
	}
	else
	{
		for (int octant = 0; octant < 8; ++octant)
		{
			const auto octant_ch = get_octant_center_halfsize(octant, my_center, my_halfsize);
			_children[octant].insert(tri_index, tri_data, octant_ch.first, octant_ch.second);
		}
	}

	return true;
}

std::vector<size_t> const & octree_child::triangles() const
{
	return _triangles;
}

int octree_child::my_octant() const
{
	assert(_parent != nullptr);
	assert(this - _parent->_children >= 0 && this - _parent->_children < 8);
	return int(this - _parent->_children);
}

void octree_child::split(vec3 my_center, vec3 my_halfsize)
{
	assert(is_leaf());
	_children = new octree_child[8];
	assert(_children != nullptr);

	std::vector<triangle> const & triangle_data = get_root().triangles();
	for (int octant = 0; octant < 8; ++octant)
	{
		_children[octant]._parent = this;
		const auto octant_ch = get_octant_center_halfsize(octant, my_center, my_halfsize);
		for (auto tri : _triangles)
			_children[octant].insert(tri, triangle_data[tri], octant_ch.first, octant_ch.second);
	}

	clear_triangles();
}

void octree_child::clear_triangles()
{
	std::vector<size_t>().swap(_triangles);
}

octree_root::octree_root(vec3 min, vec3 max)
	: center(.5*(min+max)), halfsize(.5*(max-min))
{
	_children = new octree_child[8];
	for (int octant = 0; octant < 8; ++octant)
		_children[octant]._parent = this;
}

octree_root::octree_root(std::vector<triangle> const & triangles)
{
	if (triangles.empty())
		return;

	// Find AABB_min and max
	vec3 AABB_min = triangles[0].AABB_min();
	vec3 AABB_max = triangles[0].AABB_max();

	for (const triangle t : triangles)
	{
		const vec3 tri_min = t.AABB_min();
		const vec3 tri_max = t.AABB_max();

		AABB_min =
		{
			std::min(AABB_min.x, tri_min.x),
			std::min(AABB_min.y, tri_min.y),
			std::min(AABB_min.z, tri_min.z)
		};
		AABB_max =
		{
			std::max(AABB_max.x, tri_max.x),
			std::max(AABB_max.y, tri_max.y),
			std::max(AABB_max.z, tri_max.z)
		};
	}

	// Set AABB, create child nodes
	center   = .5*(AABB_min + AABB_max);
	halfsize = .5*(AABB_max - AABB_min) + vec3{1, 1, 1};
	_children = new octree_child[8];
	for (int octant = 0; octant < 8; ++octant)
		_children[octant]._parent = this;

	// Insert triangles
	for (const auto tri : triangles)
		insert(tri);
}

int octree_root::depth() const
{
	return 0;
}

std::pair<vec3, vec3> octree_root::get_center_halfsize() const
{
	return { center, halfsize };
}

uint64_t octree_root::morton() const
{
	return 1;
}

bool octree_root::insert(triangle const & tri)
{
	if (!overlap(tri))
		return false;

	size_t tri_index = _triangles.size();
	_triangles.push_back(tri);

	for (int octant = 0; octant < 8; ++octant)
	{
		auto octant_ch = get_octant_center_halfsize(octant);
		_children[octant].insert(tri_index, tri, octant_ch.first, octant_ch.second);
	}

	return true;
}

std::vector<triangle> const & octree_root::triangles() const
{
	return _triangles;
}

linearized_octree::linearized_octree()
{
}

linearized_octree::linearized_octree(octree_root const & root) :
	center(root.center),
	halfsize(root.halfsize)
{
	if (root.triangles().size() > std::numeric_limits<int>::max())
		throw std::runtime_error("Too many triangles for linearized octree.");

	// Get a vector of all octree nodes.
	// NOTE: because of the order in which we iterate, they are sorted by Morton order.
	std::vector<octree_node const *> morton_nodes;
	{
		std::queue<octree_node const *> temp_queue;
		temp_queue.push(&root);
		while (!temp_queue.empty())
		{
			auto this_node = temp_queue.front();
			temp_queue.pop();
			morton_nodes.push_back(this_node);
			if (!this_node->is_leaf())
			{
				for (int octant = 0; octant < 8; ++octant)
					temp_queue.push(this_node->get_child(octant));
			}
		}
	}
	// Verify that Morton ordering is obeyed.
	assert(std::is_sorted(morton_nodes.begin(), morton_nodes.end(),
		[](octree_node const * a, octree_node const * b) -> bool
		{ return a->morton() < b->morton(); }));


	// Sort the triangles
	std::vector<int> triangle_map(root._triangles.size()); // For an index in root.triangles, points to the corresponding position in triangle_data.
	{
		triangle_data.reserve(root._triangles.size());

		// Temp vector, keeps track of whether a given triangle is already in triangle_data.
		std::vector<bool> done(root._triangles.size(), false);
		for (auto this_node : morton_nodes)
		{
			if (!this_node->is_leaf())
				continue;

			for (size_t tri_idx : dynamic_cast<octree_child const *>(this_node)->triangles())
			{
				if (done[tri_idx])
					continue;

				done[tri_idx] = true;
				triangle_map[tri_idx] = int(triangle_data.size());
				triangle_data.push_back(root._triangles[tri_idx]);
			}
		}

		// Verify that all triangles in root were visited.
		assert(std::find(done.begin(), done.end(), false) == done.end());
	}


	// Build linearized octree index table
	std::map<const octree_node*, int> node_p_map; // map from node pointers to their start in octree_data
	for (auto this_node : morton_nodes)
	{
		node_p_map[this_node] = int(octree_data.size());
		if (this_node->is_leaf())
		{
			// Leaf? Add triangle indices, and -1 at the end.
			// (Indices to the SORTED triangles vector, that is.)
			for(size_t tri_idx : dynamic_cast<octree_child const *>(this_node)->triangles())
				octree_data.push_back(triangle_map[tri_idx]);
			octree_data.push_back(-1);
		}
		else
		{
			// Not leaf? Add 8 empty elements for children.
			// Values will be filled in later.
			for(int octant = 0; octant < 8; octant++)
				octree_data.push_back(0);
		}
	}
	for(auto nip : node_p_map)
	{
		const octree_node* node_p = nip.first;
		const size_t index = nip.second;

		if(node_p->is_leaf())
			continue;

		for(int octant = 0; octant < 8; octant++)
		{
			const octree_child* child_p = node_p->get_child(octant);
			octree_data[index+octant] = node_p_map[child_p];
			if (child_p->is_leaf())
				octree_data[index+octant] *= -1;
		}
	}

	// If this has happened, some of the indices in octree_data (which are ints) are invalid
	if (octree_data.size() > std::numeric_limits<int>::max())
		throw std::runtime_error("Octree is too large.");
}

vec3 linearized_octree::AABB_min() const
{
	return center - halfsize;
}

vec3 linearized_octree::AABB_max() const
{
	return center + halfsize;
}

}}} // namespace nbl::geometry::octree_builder
