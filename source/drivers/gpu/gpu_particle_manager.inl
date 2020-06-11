#include <vector>
// Kill warnings for cub
#pragma warning(push, 0)
#include <cub/cub/cub.cuh>
#pragma warning(pop)

namespace nbl { namespace drivers {

template<typename material_manager_t>
CPU gpu_particle_manager<material_manager_t>
	gpu_particle_manager<material_manager_t>::create(particle_index_t capacity)
{
	gpu_particle_manager<material_manager_t> manager;
	manager._capacity = capacity;

	/*
	 * Allocate and fill particle data
	 */
	cuda::cuda_new<status_t>(&manager._status, capacity);
	cuda::cuda_new<particle_index_t>(&manager._particle_idx, capacity);
	cuda::cuda_new<particle>(&manager._particles, capacity);
	cuda::cuda_new<uint32_t>(&manager._tags, capacity);
	cuda::cuda_new<material_index_t>(&manager._material_idx, capacity);
	cuda::cuda_new<triangle*>(&manager._last_triangle, capacity);

	// Fill with initial values
	cuda::cuda_mem_scope<status_t>(manager._status, capacity,
		[capacity](status_t* device)
		{
			for (particle_index_t i = 0; i < capacity; ++i)
				device[i] = TERMINATED;
		});
	cuda::cuda_mem_scope<particle_index_t>(manager._particle_idx, capacity,
		[capacity](particle_index_t* device)
		{
			for (particle_index_t i = 0; i < capacity; ++i)
				device[i] = i;
		});
	cuda::cuda_mem_scope<material_index_t>(manager._material_idx, capacity,
		[capacity](material_index_t* device)
		{
			for (particle_index_t i = 0; i < capacity; ++i)
				device[i] = -123; // TODO: VACUUM symbolic name
		});
	cuda::cuda_mem_scope<triangle*>(manager._last_triangle, capacity,
		[capacity](triangle** device)
		{
			for (particle_index_t i = 0; i < capacity; ++i)
				device[i] = nullptr;
		});


	/*
	 * Init sorting
	 */
	// Allocate radix_index, filled with 0, 1, ... capacity
	cuda::cuda_new<particle_index_t>(&manager._radix_index, capacity);
	cuda::cuda_mem_scope<particle_index_t>(manager._radix_index, capacity,
		[&](particle_index_t* index_p)
		{
			for (particle_index_t i = 0; i < capacity; ++i)
				index_p[i] = i;
		});

	// Dump sorted status codes here
	cuda::cuda_new<status_t>(&manager._radix_dump, capacity);

	// Request radix_temp_size
	cub::DeviceRadixSort::SortPairs<status_t, particle_index_t>(
		nullptr, manager._radix_temp_size,
		manager._status, manager._radix_dump,
		manager._radix_index, manager._particle_idx,
		capacity
	);

	// Allocate radix_temp_size
	cuda::cuda_new<uint8_t>((uint8_t**)&manager._radix_temp, manager._radix_temp_size);
	return manager;
}

template<typename material_manager_t>
CPU void gpu_particle_manager<material_manager_t>::destroy(
	gpu_particle_manager<material_manager_t> & manager)
{
	cudaFree(manager._status);
	cudaFree(manager._particle_idx);
	cudaFree(manager._particles);
	cudaFree(manager._material_idx);
	cudaFree(manager._last_triangle);
	cudaFree(manager._radix_temp);
	cudaFree(manager._radix_dump);
	cudaFree(manager._radix_index);

	manager._capacity = 0;
	manager._status = nullptr;
	manager._particle_idx = nullptr;
	manager._particles = nullptr;
	manager._material_idx = nullptr;
	manager._last_triangle = nullptr;
	manager._radix_temp = nullptr;
	manager._radix_dump = nullptr;
	manager._radix_index = nullptr;
}

template<typename material_manager_t>
CPU auto gpu_particle_manager<material_manager_t>::push(
	particle* particles, uint32_t* tags, particle_index_t N)
-> particle_index_t
{
	const auto capacity = _capacity;

	// Get indices to free slots
	std::vector<particle_index_t> free_indices;
	cuda::cuda_mem_scope<status_t>(_status, _capacity,
		[&free_indices, capacity, N](status_t* status_p)
		{
			for(particle_index_t i = 0; i < capacity; ++i)
			{
				if(free_indices.size() == N)
					break;

				if(status_p[i] == TERMINATED)
				{
					free_indices.push_back(i);
				}
			}
		});

	// Copy data to GPU.
	// Do not set particle_idx.
	cuda::cuda_mem_scope<status_t>(_status, _capacity,
		[&free_indices](status_t* status_p)
		{
			for(auto idx : free_indices)
				status_p[idx] = NO_EVENT;
		});
	cuda::cuda_mem_scope<particle>(_particles, _capacity,
		[&free_indices, particles](particle* particle_p)
		{
			for(size_t i = 0; i < free_indices.size(); ++i)
				particle_p[free_indices[i]] = particles[i];
		});
	cuda::cuda_mem_scope<uint32_t>(_tags, _capacity,
		[&free_indices, tags](uint32_t* tag_p)
		{
			for(size_t i = 0; i < free_indices.size(); ++i)
				tag_p[free_indices[i]] = tags[i];
		});
	cuda::cuda_mem_scope<material_index_t>(_material_idx, _capacity,
		[&free_indices](material_index_t* material_idx_p)
		{
			for(auto idx : free_indices)
				material_idx_p[idx] = -123;
		});
	cuda::cuda_mem_scope<triangle*>(_last_triangle, _capacity,
		[&free_indices](triangle** last_triangle_p)
		{
			for(auto idx : free_indices)
				last_triangle_p[idx] = nullptr;
		});

	// free_indices.size() can never return more than the maximum value for particle_index_t
	return static_cast<particle_index_t>(free_indices.size());
}

template<typename material_manager_t>
CPU auto gpu_particle_manager<material_manager_t>::get_running_count() const -> particle_index_t
{
	particle_index_t N = 0;
	const particle_index_t capacity = _capacity;

	cuda::cuda_mem_scope<status_t>(_status, capacity,
		[&N, capacity](status_t* status)
		{
			for (particle_index_t i = 0; i < capacity; ++i)
			{
				if (status[i] != TERMINATED && status[i] != DETECTED)
				{
					++N;
				}
			}
		});

	return N;
}

template<typename material_manager_t>
CPU auto gpu_particle_manager<material_manager_t>::get_detected_count() const -> particle_index_t
{
	particle_index_t N = 0;
	const particle_index_t capacity = _capacity;

	cuda::cuda_mem_scope<status_t>(_status, capacity,
		[&N, capacity](status_t* status)
		{
			for (particle_index_t i = 0; i < capacity; ++i)
			{
				if (status[i] == DETECTED)
					++N;
			}
		});

	return N;
}

template<typename material_manager_t>
template<typename detect_function>
CPU void gpu_particle_manager<material_manager_t>::flush_detected(detect_function func)
{
	const auto capacity = _capacity;
	const auto particle_ptr = _particles;
	const auto tag_ptr = _tags;

	cuda::cuda_mem_scope<status_t>(_status, capacity,
		[&func, capacity, particle_ptr, tag_ptr](status_t* status)
		{
	cuda::cuda_mem_scope<particle>(particle_ptr, capacity,
		[&func, capacity, status, tag_ptr](particle* parts)
		{
	cuda::cuda_mem_scope<uint32_t>(tag_ptr, capacity,
		[&func, capacity, status, parts](uint32_t* tags)
		{
			for (particle_index_t i = 0; i < capacity; ++i)
			{
				if (status[i] == DETECTED)
				{
					// Do something useful with the data, then kill the particle.
					func(parts[i], tags[i]);
					status[i] = TERMINATED;
				}
			}
		});
	});
	});
}

template<typename material_manager_t>
CPU void gpu_particle_manager<material_manager_t>::sort()
{
	cub::DeviceRadixSort::SortPairs<status_t, particle_index_t>(
		_radix_temp, _radix_temp_size,
		_status, _radix_dump,
		_radix_index, _particle_idx,
		_capacity, 0, 2
	);
}

template<typename material_manager_t>
PHYSICS particle & gpu_particle_manager<material_manager_t>::operator[](particle_index_t i)
{
	return _particles[i];
}
template<typename material_manager_t>
PHYSICS particle const & gpu_particle_manager<material_manager_t>::operator[](particle_index_t i) const
{
	return _particles[i];
}

template<typename material_manager_t>
PHYSICS auto gpu_particle_manager<material_manager_t>::get_capacity() const
-> particle_index_t
{
	return _capacity;
}
template<typename material_manager_t>
PHYSICS bool gpu_particle_manager<material_manager_t>::exists(particle_index_t i) const
{
	return i < _capacity;
}

template<typename material_manager_t>
PHYSICS bool gpu_particle_manager<material_manager_t>::active(particle_index_t i) const
{
	switch(_status[i])
	{
		case TERMINATED:
		case DETECTED:
		case PENDING:
			return false;
		default:
			return true;
	}
}

template<typename material_manager_t>
PHYSICS auto gpu_particle_manager<material_manager_t>::get_particle_index(particle_index_t i) const
-> particle_index_t
{
	return _particle_idx[i];
}

template<typename material_manager_t>
PHYSICS auto gpu_particle_manager<material_manager_t>::get_material_index(particle_index_t i) const
-> material_index_t
{
	return _material_idx[i];
}
template<typename material_manager_t>
PHYSICS void gpu_particle_manager<material_manager_t>::set_material_index(
	particle_index_t particle_idx, material_index_t new_material_idx)
{
	_material_idx[particle_idx] = new_material_idx;
}

template<typename material_manager_t>
PHYSICS triangle const * gpu_particle_manager<material_manager_t>::get_last_triangle(particle_index_t i) const
{
	return _last_triangle[i];
}
template<typename material_manager_t>
PHYSICS void gpu_particle_manager<material_manager_t>::forget_last_triangle(particle_index_t i)
{
	_last_triangle[i] = nullptr;
}

template<typename material_manager_t>
PHYSICS bool gpu_particle_manager<material_manager_t>::next_elastic(particle_index_t i) const
{
	return _status[i] == ELASTIC_EVENT;
}
template<typename material_manager_t>
PHYSICS bool gpu_particle_manager<material_manager_t>::next_inelastic(particle_index_t i) const
{
	return _status[i] == PENDING || _status[i] == INELASTIC_EVENT;
}
template<typename material_manager_t>
PHYSICS bool gpu_particle_manager<material_manager_t>::next_intersect(particle_index_t i) const
{
	return _status[i] == INTERSECT_EVENT;
}
template<typename material_manager_t>
PHYSICS bool gpu_particle_manager<material_manager_t>::is_detected(particle_index_t i) const
{
	return _status[i] == DETECTED;
}
template<typename material_manager_t>
PHYSICS bool gpu_particle_manager<material_manager_t>::is_terminated(particle_index_t i) const
{
	return _status[i] == TERMINATED;
}

template<typename material_manager_t>
PHYSICS void gpu_particle_manager<material_manager_t>::add_particle(
	particle_index_t target_idx,
	particle new_particle,
	uint32_t new_tag)
{
	_status[target_idx] = NO_EVENT;
	_particles[target_idx] = new_particle;
	_tags[target_idx] = new_tag;
	_material_idx[target_idx] = -123; // TODO
	_last_triangle[target_idx] = nullptr;
}
template<typename material_manager_t>
PHYSICS bool gpu_particle_manager<material_manager_t>::secondary_slot_free() const
{
	return _status[get_secondary_slot()] == TERMINATED;
}
template<typename material_manager_t>
PHYSICS void gpu_particle_manager<material_manager_t>::create_secondary(
	particle_index_t primary_idx, particle secondary_particle)
{
	const particle_index_t secondary_idx = get_secondary_slot();
	_status[secondary_idx] = NEW_SECONDARY;
	_particles[secondary_idx] = secondary_particle;
	_tags[secondary_idx] = _tags[primary_idx];
	_material_idx[secondary_idx] = _material_idx[primary_idx];
	_last_triangle[secondary_idx] = nullptr;
}
template<typename material_manager_t>
PHYSICS auto gpu_particle_manager<material_manager_t>::get_secondary_slot() const
-> particle_index_t
{
	particle_index_t thread_idx = threadIdx.x + blockIdx.x*blockDim.x;
	return _particle_idx[_capacity - 1 - thread_idx];
}

template<typename material_manager_t>
PHYSICS void gpu_particle_manager<material_manager_t>::terminate(particle_index_t i)
{
	_status[i] = TERMINATED;
}
template<typename material_manager_t>
PHYSICS void gpu_particle_manager<material_manager_t>::detect(
	particle_index_t i)
{
	_status[i] = DETECTED;
}

template<typename material_manager_t>
PHYSICS void gpu_particle_manager<material_manager_t>::set_scatter_event(
	particle_index_t i, uint8_t event)
{
	// TODO: event types hardcoded here
	switch (event)
	{
	case 1:
		_status[i] = INELASTIC_EVENT;
		break;
	case 2:
		_status[i] = ELASTIC_EVENT;
		break;
	default:
		_status[i] = NO_EVENT;
		break;
	}
}
template<typename material_manager_t>
PHYSICS void gpu_particle_manager<material_manager_t>::set_intersect_event(
	particle_index_t i, triangle* t)
{
	_status[i] = INTERSECT_EVENT;
	_last_triangle[i] = t;
}

template<typename material_manager_t>
PHYSICS void gpu_particle_manager<material_manager_t>::pending(particle_index_t i)
{
	_status[i] = PENDING;
}
template<typename material_manager_t>
PHYSICS void gpu_particle_manager<material_manager_t>::inelastic(particle_index_t i)
{
	_status[i] = INELASTIC_EVENT;
}

}} // namespace nbl::drivers
