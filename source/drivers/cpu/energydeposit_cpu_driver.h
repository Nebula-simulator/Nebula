#ifndef __ENERGYDEPOSIT_CPU_DRIVER_H_
#define __ENERGYDEPOSIT_CPU_DRIVER_H_

#include "cpu_driver.h"

namespace nbl { namespace drivers {

/**
 * \brief CPU driver for tracking energy deposits during scattering.
 * \see cpu_driver
 */

template<
	typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
class energydeposit_cpu_driver
	: public cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>
{
public:
	using base_t = cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>;

	using typename base_t::particle_index_t;
	using typename base_t::material_t;
	using typename base_t::material_manager_t;
	using typename base_t::seed_t;

	/// Constructor
	energydeposit_cpu_driver(
		intersect_t intersect,
		material_manager_t const & materials,
		geometry_manager_t const & geometry,
		real min_energy, real max_energy,
		seed_t seed = util::random_generator<false>::default_seed)
	: base_t(intersect, materials, geometry, min_energy, max_energy, seed)
	{}

	/**
	 * \brief Perform a single iteration of the simulation for all particles.
	 *
	 * \param function Callback function to be called if a particle loses energy
	 *                 in a scattering event. Should have signature
	 *                 `void(particle const & before, particle const & after, uint32_t tag)`.
	 */
	template<typename deposit_function>
	void do_iteration(deposit_function function)
	{
		const auto particle_count = this->_particles.get_total_count();
		for (particle_index_t particle_idx = 0; particle_idx < particle_count; ++particle_idx)
		{
			if (!this->_particles.active(particle_idx))
				continue;

			this->init(particle_idx);

			// We want to track only scattering events, not interface events
			bool track = this->_particles.next_scatter(particle_idx);
			const particle before = this->_particles[particle_idx];

			this->intersect(particle_idx);
			this->scatter(particle_idx);

			const particle after = this->_particles[particle_idx];
			if (track && before.kin_energy != after.kin_energy)
				function(before, after,
					this->_particles.get_primary_tag(particle_idx));
		}
	}

	/// Keep simulating until there are no particles left.
	template<typename deposit_function>
	void simulate_to_end(deposit_function&& func)
	{
		for (particle_index_t particle_idx = 0;
			particle_idx < this->_particles.get_total_count(); ++particle_idx)
		{
			while (this->_particles.active(particle_idx))
			{
				this->init(particle_idx);

				// We want to track only scattering events, not interface events
				bool track = this->_particles.next_scatter(particle_idx);
				const particle before = this->_particles[particle_idx];

				this->intersect(particle_idx);
				this->scatter(particle_idx);

				const particle after = this->_particles[particle_idx];
				if (track && before.kin_energy != after.kin_energy)
					func(before, after,
						this->_particles.get_primary_tag(particle_idx));
			}
		}
	}
};

}} // namespace nbl::drivers

#endif // __ENERGYDEPOSIT_CPU_DRIVER_H_
