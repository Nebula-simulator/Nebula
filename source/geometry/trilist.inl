namespace nbl { namespace geometry {

template<bool gpu_flag>
CPU trilist<gpu_flag> trilist<gpu_flag>::create(std::vector<triangle> const & triangles)
{
	// TODO: error message
	if (triangles.empty())
		throw std::runtime_error("No triangles provided!");

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

	AABB_min -= vec3{ 1, 1, 1 };
	AABB_max += vec3{ 1, 1, 1 };

	return detail::trilist_factory<gpu_flag>::create(triangles, AABB_min, AABB_max);
}

template<bool gpu_flag>
CPU void trilist<gpu_flag>::destroy(trilist<gpu_flag> & geometry)
{
	detail::trilist_factory<gpu_flag>::destroy(geometry);
}

template<bool gpu_flag>
PHYSICS bool trilist<gpu_flag>::in_domain(vec3 pos)
{
	return ((pos.x > _AABB_min.x) && (pos.x < _AABB_max.x)
		&& (pos.y > _AABB_min.y) && (pos.y < _AABB_max.y)
		&& (pos.z > _AABB_min.z) && (pos.z < _AABB_max.z));
}

template<bool gpu_flag>
PHYSICS intersect_event trilist<gpu_flag>::propagate(vec3 start, vec3 direction, real distance,
	triangle const * ignore_triangle, int ignore_material) const
{
	intersect_event evt { distance, nullptr };
	for (triangle_index_t i = 0; i < _N; ++i)
	{
		if (_triangles + i == ignore_triangle)
			continue;

		// retrieve triangle from global memory
		const triangle this_triangle = _triangles[i];

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

		const real t = this_triangle.intersect_ray(start, direction);
		if (t > 0 && t < evt.isect_distance)
		{
			evt.isect_distance = t;
			evt.isect_triangle = _triangles + i;
		}
	}
	return evt;
}

template<bool gpu_flag>
PHYSICS real trilist<gpu_flag>::get_max_extent() const
{
	return _max_extent;
}

template<bool gpu_flag>
inline PHYSICS vec3 trilist<gpu_flag>::AABB_min() const
{
	return _AABB_min;
}
template<bool gpu_flag>
inline PHYSICS vec3 trilist<gpu_flag>::AABB_max() const
{
	return _AABB_max;
}

template<bool gpu_flag>
CPU void trilist<gpu_flag>::set_AABB(vec3 min, vec3 max)
{
	_AABB_min = min;
	_AABB_max = max;

	const vec3 m = max - min;
	_max_extent = magnitude(m);
}


namespace detail
{
	template<>
	struct trilist_factory<false>
	{
		inline static CPU trilist<false> create(std::vector<triangle> triangles, vec3 AABB_min, vec3 AABB_max)
		{
			using trilist_t = trilist<false>;
			using triangle_index_t = trilist_t::triangle_index_t;

			trilist_t geometry;

			if (triangles.size() > std::numeric_limits<triangle_index_t>::max())
				throw std::runtime_error("Too many triangles in geometry");
			geometry._N = static_cast<triangle_index_t>(triangles.size());

			geometry._triangles = reinterpret_cast<triangle*>(malloc(geometry._N * sizeof(triangle)));
			for (triangle_index_t i = 0; i < triangles.size(); ++i)
			{
				geometry._triangles[i] = triangles[i];
			}

			geometry.set_AABB(AABB_min, AABB_max);

			return geometry;
		}

		inline static CPU void free(trilist<false> & geometry)
		{
			::free(geometry._triangles);

			geometry._triangles = nullptr;
			geometry._N = 0;
		}
	};

#if CUDA_COMPILER_AVAILABLE
	template<>
	struct trilist_factory<true>
	{
		inline static CPU trilist<true> create(std::vector<triangle> triangles, vec3 AABB_min, vec3 AABB_max)
		{
			using trilist_t = trilist<true>;
			using triangle_index_t = trilist_t::triangle_index_t;

			trilist_t geometry;

			if (triangles.size() > std::numeric_limits<triangle_index_t>::max())
				throw std::runtime_error("Too many triangles in geometry");
			geometry._N = static_cast<triangle_index_t>(triangles.size());

			// Copy triangle data to device
			cuda::cuda_new<triangle>(&geometry._triangles, geometry._N);
			cuda::cuda_mem_scope<triangle>(geometry._triangles, geometry._N,
				[&triangles](triangle* device)
			{
				for (triangle_index_t i = 0; i < triangles.size(); ++i)
					device[i] = triangles[i];
			});

			geometry.set_AABB(AABB_min, AABB_max);

			return geometry;
		}

		inline static CPU void free(trilist<true> & geometry)
		{
			cudaFree(geometry._triangles);

			geometry._triangles = nullptr;
			geometry._N = 0;
		}
	};
#endif // CUDA_COMPILER_AVAILABLE
} // namespace detail

}} // namespace nbl::geometry
