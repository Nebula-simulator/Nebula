#ifndef __CPU_DRIVER_H_
#define __CPU_DRIVER_H_

#include "../../core/cpu_material_manager.h"
#include "../../common/util/random.h"
#include "cpu_particle_manager.h"

namespace nbl { namespace drivers {

/**
 * \brief Generic CPU driver
 *
 * Accepts any combination of scattering mechanisms, any number of secondaries
 * may be generated per event.
 *
 * This is a base class without any functions to perform a simulation. "Useful"
 * drivers inherit from this class and provide such functionality.
 * Such drivers may, for instance, keep track of scattering events.
 *
 * \tparam scatter_list_t     List of scattering mechanisms. Should be of type
 *                            ::scatter_list.
 * \tparam intersect_t        Intersection event handler
 * \tparam geometry_manager_t Geometry manager
 */
template<
	typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
class cpu_driver
{
public:
	using material_t = material<scatter_list_t>;
	using material_manager_t = cpu_material_manager<material_t>;
	using particle_manager_t = cpu_particle_manager<material_manager_t>;

	using particle_index_t = typename particle_manager_t::particle_index_t;
	using primary_tag_t = typename particle_manager_t::primary_tag_t;
	using seed_t = typename util::random_generator<false>::seed_type;

	/**
	 * \brief Add new primary particles to the simulation.
	 *
	 * \param particles Pointer to an array of particles to be added.
	 * \param tags      Pointer to the corresponding array of tags.
	 * \param N         Number of particles to be added.
	 */
	inline particle_index_t push(particle* particles, primary_tag_t* tags, particle_index_t N);

	/**
	 * \brief Get number of particles currently in the simulation.
	 *
	 * More specifically, those that are not terminated or detected.
	 */
	particle_index_t get_running_count() const;

	/// Get number of detected particles currently in the simulation.
	particle_index_t get_detected_count() const;

	/**
	 * \brief Remove detected particles from the simulation, calling a callback
	 *        before doing so.
	 *
	 * The callback function receives the particle data and the tag that
	 * belonged to the primary electron that initated the cascade the detected
	 * particle is part of.
	 *
	 * \param function Callback function to be called for each detected particle.
	 *                 Should have signature void(particle const &, uint32_t).
	 */
	template<typename detect_function>
	void flush_detected(detect_function function);

	// Energy threshold (w.r.t. the vacuum level) below which particles must be terminated.
	real energy_threshold;

protected:
	/**
	 * \brief Constructor.
	 *
	 * \param intersect        Instance of the intersection handler.
	 * \param materials        List of materials in the simulation.
	 * \param geometry         Geometry manager, holding the simulation geometry.
	 * \param energy_threshold Energy threshold w.r.t. the vacuum level below
	 *                         which particles must be terminated.
	 * \param seed             Seed for the random number generator.
	 */
	cpu_driver(
		intersect_t intersect,
		material_manager_t const & materials,
		geometry_manager_t const & geometry,
		real energy_threshold = 0,
		seed_t seed = util::random_generator<false>::default_seed);

	/// Destructor
	~cpu_driver();

	particle_manager_t _particles;
	material_manager_t const & _materials;
	geometry_manager_t const & _geometry;
	intersect_t _intersect;

	// Random number generator
	util::random_generator<false> rand_state;

	// Functions used in simulations
	inline void init(particle_index_t particle_idx);
	inline void intersect(particle_index_t particle_idx);
	inline void scatter(particle_index_t particle_idx);

	cpu_driver(cpu_driver const &) = delete;
	cpu_driver& operator=(cpu_driver const &) = delete;
};

}} // namespace nbl::drivers

#include "cpu_driver.inl"

#endif // __CPU_DRIVER_H_
