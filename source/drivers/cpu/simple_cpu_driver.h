#ifndef __SIMPLE_CPU_DRIVER_H_
#define __SIMPLE_CPU_DRIVER_H_

#include "cpu_driver.h"

namespace nbl { namespace drivers {

/**
 * \brief CPU driver that only tries to perform the simulation and nothing else.
 * \see cpu_driver
 */

template<
	typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
class simple_cpu_driver
	: public cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>
{
public:
	using base_t = cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>;

	using typename base_t::particle_index_t;
	using typename base_t::material_t;
	using typename base_t::seed_t;

	/// Constructor
	simple_cpu_driver(
		geometry_manager_t const & geometry,
		intersect_t const & intersect,
		std::vector<material_t> const & materials,
		real energy_threshold = 0,
		seed_t seed = util::random_generator<false>::default_seed)
	: base_t(geometry, intersect, materials, energy_threshold, seed)
	{}

	/// Perform a single iteration of the simulation for all particles.
	void do_iteration()
	{
		const auto particle_count = this->_particles.get_total_count();
		for (particle_index_t particle_idx = 0; particle_idx < particle_count; ++particle_idx)
		{
			if (!this->_particles.active(particle_idx))
				continue;

			this->init(particle_idx);
			this->intersect(particle_idx);
			this->scatter(particle_idx);
		}
	}

	/// Keep simulating until there are no particles left.
	void simulate_to_end()
	{
		for (particle_index_t particle_idx = 0;
			particle_idx < this->_particles.get_total_count(); ++particle_idx)
		{
			while (this->_particles.active(particle_idx))
			{
				this->init(particle_idx);
				this->intersect(particle_idx);
				this->scatter(particle_idx);
			}
		}
	}
};

}} // namespace nbl::drivers

#endif // __SIMPLE_CPU_DRIVER_H_
