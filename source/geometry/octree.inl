namespace nbl { namespace geometry {

template<bool gpu_flag>
CPU octree<gpu_flag> octree<gpu_flag>::create(octree_builder::linearized_octree const & linoct)
{
	return detail::octree_factory<gpu_flag>::create(linoct);
}

template<bool gpu_flag>
CPU octree<gpu_flag> octree<gpu_flag>::create(std::vector<triangle> const & triangles)
{
	return detail::octree_factory<gpu_flag>::create(
		octree_builder::linearized_octree(
			octree_builder::octree_root(triangles)
		));
}

template<bool gpu_flag>
CPU void octree<gpu_flag>::destroy(octree<gpu_flag> & geometry)
{
	detail::octree_factory<gpu_flag>::free(geometry);
}

template<bool gpu_flag>
PHYSICS bool octree<gpu_flag>::in_domain(vec3 pos) const
{
	return ((pos.x > _AABB_center.x-_AABB_halfsize.x) && (pos.x < _AABB_center.x+_AABB_halfsize.x)
		&& (pos.y > _AABB_center.y-_AABB_halfsize.y) && (pos.y < _AABB_center.y+_AABB_halfsize.y)
		&& (pos.z > _AABB_center.z-_AABB_halfsize.z) && (pos.z < _AABB_center.z+_AABB_halfsize.z));
}

template<bool gpu_flag>
PHYSICS intersect_event octree<gpu_flag>::propagate(vec3 start, vec3 direction, real distance,
	triangle const * ignore_triangle, int ignore_material) const
{
	// Current (x, y, z) position in the octree, is updated as we travel through cells
	vec3 current_position = start;
	// Distance we still need to travel from current_position, in units of direction
	real distance_to_go = distance;

	// Search in octree
	// A cell in the octree is uniquely identified by a location code.
	// The root cell has location code 1.
	// The location code of a decedant is obtained by shifting the location code
	// to the left by three bits and adding the octant of the decadant to the location code.
	uint64_t location = 1;
	do {
		// define the root axis aligned bounding box
		vec4 AABB{
			_AABB_center.x,
			_AABB_center.y,
			_AABB_center.z,
			1 // scale factor
		};

		// initial index to linearized octree
		int index = 0;

		// traverse to location: we compute the location of the current node and
		//   set the index to the node in the linearized octree array
		// `__clzll` computes number of leading zeros, effectively this is 64-floor(log2(n))
		for (int i = 60 - clz(location); i >= 0; i -= 3) {
			const unsigned int octant = (location >> i) & 7;
			index = _octree_data[index + octant];
			AABB.x += AABB.w*_AABB_halfsize.x*((octant & 1) - 0.5_r);
			AABB.y += AABB.w*_AABB_halfsize.y*(((octant & 2) >> 1) - 0.5_r);
			AABB.z += AABB.w*_AABB_halfsize.z*(((octant & 4) >> 2) - 0.5_r);
			AABB.w *= 0.5_r;
		}

		// traverse to leaf node (which has a linearized octree index smaller than zero)
		while (index >= 0) {
			unsigned int octant = 0;
			octant ^= (current_position.x > AABB.x) ? 1 : 0;
			octant ^= (current_position.y > AABB.y) ? 2 : 0;
			octant ^= (current_position.z > AABB.z) ? 4 : 0;
			location = (location << 3) | octant;
			index = _octree_data[index + octant];
			AABB.x += AABB.w*_AABB_halfsize.x*((octant & 1) - 0.5_r);
			AABB.y += AABB.w*_AABB_halfsize.y*(((octant & 2) >> 1) - 0.5_r);
			AABB.z += AABB.w*_AABB_halfsize.z*(((octant & 4) >> 2) - 0.5_r);
			AABB.w *= 0.5_r;
		}
		// note that a leaf has a guaranteed negative index, the actual index to the leaf
		// node is then obtained by negating.
		index = -index;


		// determine cell/triangle intersections for the leaf
		real intersect; // Distance to next triangle intersection
		int target; // target >= 0 : triangle index
					// target = -1 : cell x intersection
					// target = -2 : cell y intersection
					// target = -4 : cell z intersection

		// find intersection times with each orthogonal plane of the leaf's bounding box,
		// then see which plane is reached first.
		{
			const vec3 t = AABB_intersect(current_position, direction, { AABB.x, AABB.y, AABB.z }, AABB.w * _AABB_halfsize);
			if ((t.x < t.y) && (t.x < t.z)) {
				intersect = t.x;
				target = -1;
			}
			else if ((t.y < t.x) && (t.y < t.z)) {
				intersect = t.y;
				target = -2;
			}
			else {
				intersect = t.z;
				target = -4;
			}
		}


		// iterate triangles in leaf
		while (true) {
			const int triangle_idx = _octree_data[index++];
			if (triangle_idx < 0)
				break;

			// don't intersect with ignore_triangle
			if (_triangles + triangle_idx == ignore_triangle)
				continue;

			// retrieve triangle from global memory
			const triangle this_triangle = _triangles[triangle_idx];

			int mat_idx_out;
			if (dot_product(this_triangle.get_normal(), direction) < 0)
				mat_idx_out = this_triangle.material_in;
			else
				mat_idx_out = this_triangle.material_out;

			// if the outgoing material is the same as current, nothing happens
			// if the triangle represents a detector which can't see the current
			// particle, nothing happens
			if (mat_idx_out == ignore_material)
				continue;

			// compute the intersection with the triangle; keep it if it
			// is closer than the current distance to next cell or triangle
			const real t = this_triangle.intersect_ray(current_position, direction);
			if ((t > 0) && (t <= intersect + EPSILON)) {
				intersect = t;
				target = triangle_idx;
			}
		}


		// manage calculated intersections
		if (intersect >= distance_to_go) {
			// EXIT: distance traveled without triangle intersections
			return { distance, nullptr };
		}
		else if (target >= 0) {
			// EXIT: triangle intersection
			return { distance - distance_to_go + intersect, _triangles + target };
		}

		// go to edge of the leaf
		distance_to_go -= intersect;
		current_position += direction * intersect;

		// find adjacent node by bit-magic, and loop
		// the adjacent node is determined solely by location code!
		unsigned int mask = -target;
		unsigned int value;
		// mask = 1 (001), value = 0 (000) : xx0 --> xx1 (find neighbor in positive x direction)
		// mask = 1 (001), value = 1 (000) : xx1 --> xx0 (find neighbor in negative x direction)
		// mask = 2 (010), value = 0 (000) : x0x --> x1x (find neighbor in positive y direction)
		// mask = 2 (010), value = 2 (010) : x1x --> x0x (find neighbor in negative y direction)
		// mask = 4 (100), value = 0 (000) : 0xx --> 1xx (find neighbor in positive z direction)
		// mask = 4 (100), value = 4 (100) : 1xx --> 0xx (find neighbor in negative z direction)
		if (mask == 1)
			value = (direction.x >= 0) ? 0 : 1;
		else if (mask == 2)
			value = (direction.y >= 0) ? 0 : 2;
		else
			value = (direction.z >= 0) ? 0 : 4;
		while (location > 1) {
			if ((location&mask) == value) {
				location ^= mask;
				break;
			}
			location >>= 3;
		}

	} while (location > 1);

	// EXIT: particle is out of grid (left root cell)
	return { distance, nullptr };
}

template<bool gpu_flag>
PHYSICS real octree<gpu_flag>::get_max_extent() const
{
	return _max_extent;
}

template<bool gpu_flag>
inline PHYSICS vec3 octree<gpu_flag>::AABB_min() const
{
	return _AABB_center - _AABB_halfsize;
}
template<bool gpu_flag>
inline PHYSICS vec3 octree<gpu_flag>::AABB_max() const
{
	return _AABB_center + _AABB_halfsize;
}

template<bool gpu_flag>
PHYSICS vec3 octree<gpu_flag>::AABB_intersect(vec3 pos, vec3 dir, vec3 center, vec3 halfsize)
{
	return
	{
		(center.x + copysignr(halfsize.x + EPSILON, dir.x) - pos.x) / dir.x,
		(center.y + copysignr(halfsize.y + EPSILON, dir.y) - pos.y) / dir.y,
		(center.z + copysignr(halfsize.z + EPSILON, dir.z) - pos.z) / dir.z
	};
}

#ifdef _MSC_VER
#include <intrin.h>
#endif

template<bool gpu_flag>
PHYSICS int octree<gpu_flag>::clz(uint64_t x)
{
	static_assert(sizeof(x) == sizeof(long long int),
		"__clzll intrinsic should be 64 bit?");
#if defined __CUDA_ARCH__
	return __clzll(x);
#elif defined __GNUC__
	return __builtin_clzll(x);
#elif defined(_MSC_VER)
	return (int)__lzcnt64(x);
#else
	static_assert(0, "No clz implementation for your platform!");
#endif
}

namespace detail
{
	template<>
	struct octree_factory<false>
	{
		inline static CPU octree<false> create(
			octree_builder::linearized_octree const & linoct)
		{
			octree<false> geometry;

			// Copy octree data
			geometry._octree_data = new int[linoct.octree_data.size()];
			memcpy(geometry._octree_data, linoct.octree_data.data(), linoct.octree_data.size() * sizeof(int));

			// Copy triangle data
			geometry._triangles = reinterpret_cast<triangle*>(malloc(linoct.triangle_data.size() * sizeof(triangle)));
			memcpy(geometry._triangles, linoct.triangle_data.data(), linoct.triangle_data.size() * sizeof(triangle));

			// Copy AABB
			geometry._AABB_center = linoct.center;
			geometry._AABB_halfsize = linoct.halfsize;
			geometry._max_extent = 2*magnitude(linoct.halfsize);

			return geometry;
		}

		inline static CPU void free(octree<false> & geometry)
		{
			delete[] geometry._octree_data;
			::free(geometry._triangles);

			geometry._octree_data = nullptr;
			geometry._triangles = nullptr;
		}
	};

#if CUDA_COMPILER_AVAILABLE
	template<>
	struct octree_factory<true>
	{
		inline static CPU octree<true> create(
			octree_builder::linearized_octree const & linoct)
		{
			octree<true> geometry;

			// Copy octree data
			cuda::cuda_new<int>(&geometry._octree_data, linoct.octree_data.size());
			cudaMemcpy(geometry._octree_data, linoct.octree_data.data(),
				linoct.octree_data.size()*sizeof(int), cudaMemcpyHostToDevice);

			// Copy triangle data
			cuda::cuda_new<triangle>(&geometry._triangles, linoct.triangle_data.size());
			cudaMemcpy(geometry._triangles, linoct.triangle_data.data(),
				linoct.triangle_data.size()*sizeof(triangle), cudaMemcpyHostToDevice);

			// Copy AABB
			geometry._AABB_center = linoct.center;
			geometry._AABB_halfsize = linoct.halfsize;
			geometry._max_extent = 2*magnitude(linoct.halfsize);

			return geometry;
		}

		inline static CPU void free(octree<true> & geometry)
		{
			cudaFree(geometry._octree_data);
			cudaFree(geometry._triangles);

			geometry._octree_data = nullptr;
			geometry._triangles = nullptr;
		}
	};
#endif // CUDA_COMPILER_AVAILABLE
} // namespace detail

}} // namespace nbl::geometry
