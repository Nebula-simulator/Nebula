#include <algorithm>

namespace nbl { namespace drivers {

template<typename material_manager_t>
cpu_particle_manager<material_manager_t>
	cpu_particle_manager<material_manager_t>::create()
{
	cpu_particle_manager<material_manager_t> manager{};
	return manager;
}

template<typename material_manager_t>
void cpu_particle_manager<material_manager_t>::destroy(
	cpu_particle_manager<material_manager_t> & manager)
{
	manager.particles.clear();
	manager.particles.shrink_to_fit();
	manager.cascades.clear();
}

template<typename material_manager_t>
auto cpu_particle_manager<material_manager_t>::push(
	particle* primary_particles, primary_tag_t* tags, particle_index_t N)
-> particle_index_t
{
	particles.reserve(particles.size() + N);
	for (particle_index_t i = 0; i < N; ++i)
	{
		particles.push_back({
			NO_EVENT,
			0,
			-123,  // TODO: vacuum
			primary_particles[i],
			tags[i],
			0,
			nullptr
		});

		cascades.insert(std::make_pair(tags[i], cascade_struct{1, 1}));
	}
	return N;
}

template<typename material_manager_t>
template<typename detect_function>
void cpu_particle_manager<material_manager_t>::flush_detected(detect_function func)
{
	for (auto& this_particle : particles)
	{
		if (this_particle.status == DETECTED)
		{
			func(this_particle.particle_data, this_particle.primary_tag);
			this_particle.status = TERMINATED;
		}
	}
}

template<typename material_manager_t>
void cpu_particle_manager<material_manager_t>::flush_terminated()
{
	particles.erase
	(
		std::remove_if(particles.begin(), particles.end(),
			[](particle_struct const & x) -> bool
			{ return x.status == TERMINATED; }),
		particles.end()
	);

	// C++20: can use std::erase_if here.
	for (auto i = cascades.begin(); i != cascades.end(); )
	{
		if (i->second.running_count == 0)
			i = cascades.erase(i);
		else
			++i;
	}
}

template<typename material_manager_t>
auto cpu_particle_manager<material_manager_t>::get_total_count() const -> particle_index_t
{
	return particles.size();
}
template<typename material_manager_t>
auto cpu_particle_manager<material_manager_t>::get_running_count() const -> particle_index_t
{
	return static_cast<particle_index_t>(
	std::count_if(particles.begin(), particles.end(),
		[](particle_struct const & x) -> bool
		{ return x.status != TERMINATED && x.status != DETECTED; }));
}
template<typename material_manager_t>
auto cpu_particle_manager<material_manager_t>::get_detected_count() const -> particle_index_t
{
	return std::count_if(particles.begin(), particles.end(),
		[](particle_struct const & x) -> bool
		{ return x.status == DETECTED; });
}

template<typename material_manager_t>
PHYSICS particle & cpu_particle_manager<material_manager_t>::operator[](particle_index_t i)
{
	return particles[i].particle_data;
}
template<typename material_manager_t>
PHYSICS particle const & cpu_particle_manager<material_manager_t>::operator[](particle_index_t i) const
{
	return particles[i].particle_data;
}

template<typename material_manager_t>
PHYSICS bool cpu_particle_manager<material_manager_t>::exists(particle_index_t i) const
{
	return i < particles.size();
}

template<typename material_manager_t>
PHYSICS bool cpu_particle_manager<material_manager_t>::active(
	particle_index_t i) const
{
	switch (particles[i].status)
	{
		case TERMINATED:
		case DETECTED:
			return false;
		default:
			return true;
	}
}

template<typename material_manager_t>
PHYSICS auto cpu_particle_manager<material_manager_t>::get_material_index(particle_index_t i) const
-> material_index_t
{
	return particles[i].current_material;
}
template<typename material_manager_t>
PHYSICS void cpu_particle_manager<material_manager_t>::set_material_index(
	particle_index_t particle_idx, material_index_t new_material_idx)
{
	particles[particle_idx].current_material = new_material_idx;
}

template<typename material_manager_t>
PHYSICS auto cpu_particle_manager<material_manager_t>::get_primary_tag(particle_index_t i) const
-> primary_tag_t
{
	return particles[i].primary_tag;
}

template<typename material_manager_t>
PHYSICS triangle const * cpu_particle_manager<material_manager_t>::get_last_triangle(
	particle_index_t i) const
{
	return particles[i].last_triangle;
}
template<typename material_manager_t>
PHYSICS void cpu_particle_manager<material_manager_t>::forget_last_triangle(
	particle_index_t i)
{
	particles[i].last_triangle = nullptr;
}

template<typename material_manager_t>
PHYSICS bool cpu_particle_manager<material_manager_t>::next_scatter(particle_index_t i) const
{
	return particles[i].status == SCATTER_EVENT;
}
template<typename material_manager_t>
PHYSICS uint8_t cpu_particle_manager<material_manager_t>::get_next_scatter(particle_index_t i) const
{
	return particles[i].next_scatter;
}
template<typename material_manager_t>
PHYSICS bool cpu_particle_manager<material_manager_t>::next_intersect(particle_index_t i) const
{
	return particles[i].status == INTERSECT_EVENT;
}

template<typename material_manager_t>
PHYSICS void cpu_particle_manager<material_manager_t>::create_secondary(
	particle_index_t primary_idx, particle secondary_particle)
{
	const auto primary_tag = particles[primary_idx].primary_tag;
	particles.push_back({
		NO_EVENT,
		0,
		get_material_index(primary_idx),
		secondary_particle,
		primary_tag,
		cascades[primary_tag].next_secondary_tag++,
		nullptr
	});
	++cascades[primary_tag].running_count;
}

template<typename material_manager_t>
PHYSICS void cpu_particle_manager<material_manager_t>::terminate(particle_index_t i)
{
	particles[i].status = TERMINATED;
	--cascades[particles[i].primary_tag].running_count;
}
template<typename material_manager_t>
PHYSICS void cpu_particle_manager<material_manager_t>::detect(particle_index_t i)
{
	particles[i].status = DETECTED;
	--cascades[particles[i].primary_tag].running_count;
}

// TODO: the two functions below recalculate the normalization of "dir"
// which has already been done by the driver...
template<typename material_manager_t>
PHYSICS void cpu_particle_manager<material_manager_t>::set_scatter_event(
	particle_index_t i, scatter_event event)
{
	if (event.type != 0)
	{
		particles[i].status = SCATTER_EVENT;
		particles[i].next_scatter = event.type;
	}
	else
	{
		particles[i].status = NO_EVENT;
	}
	particles[i].particle_data.pos += normalised(particles[i].particle_data.dir) * event.distance;
}
template<typename material_manager_t>
PHYSICS void cpu_particle_manager<material_manager_t>::set_intersect_event(
	particle_index_t i, intersect_event event)
{
	particles[i].status = INTERSECT_EVENT;
	particles[i].last_triangle = event.isect_triangle;
	particles[i].particle_data.pos += normalised(particles[i].particle_data.dir) * event.isect_distance;
}

}} // namespace nbl::drivers
