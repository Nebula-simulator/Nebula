#ifndef __GPU_DRIVER_H_
#define __GPU_DRIVER_H_

#include "../../core/material_manager.h"
#include "../../common/util/random.h"
#include "../../common/work_pool.h"
#include "gpu_particle_manager.h"

namespace nbl { namespace drivers {

/**
 * \brief Driver similar to Thomas Verduin's e-scatter.
 *
 * It only accepts two scattering mechanisms. The first of these may generate
 * secondary electrons, the second may not. Only one SE may be generated per
 * event.
 *
 * T.V.'s thesis: doi:10.4233/uuid:f214f594-a21f-4318-9f29-9776d60ab06c
 *
 * \tparam scatter_list_t     List of scattering mechanisms. Should be of type
 *                            ::scatter_list.
 * \tparam intersect_t        Intersection event handler
 * \tparam geometry_manager_t Geometry manager
 */
template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t
>
class gpu_driver
{
public:
	static_assert(scatter_list_t::size() == 2,
		"GPU driver must have two scattering mechanisms");
	static_assert(scatter_list_t::template type_at_index<1>::may_create_se == false,
		"Only the first scattering mechanism may create secondaries in GPU code");

	using material_t = material<scatter_list_t>;
	using material_manager_t = material_manager<material_t, true>;
	using particle_manager_t = gpu_particle_manager<material_manager_t>;

	using particle_index_t = typename particle_manager_t::particle_index_t;
	using status_t = typename particle_manager_t::status_t;
	using seed_t = typename util::random_generator<true>::seed_type;

	/**
	 * \brief Constructor.
	 *
	 * \param particle_capacity Maximum number of particles the simulation may have.
	 * \param geom              Geometry manager, holding the simulation geometry.
	 * \param inter             Instance of the intersection handler.
	 * \param materials         List of materials in the simulation.
	 * \param energy_threshold  Energy threshold w.r.t. the vacuum level below
	 *                          which particles must be terminated.
	 * \param seed              Seed for the random number generator.
	 */
	CPU gpu_driver(particle_index_t particle_capacity, geometry_manager_t geom,
		intersect_t inter, std::vector<material_t> materials,
		real energy_threshold = 0,
		seed_t seed = util::random_generator<true>::default_seed);
	CPU ~gpu_driver();

	/**
	 * \brief Add new primary particles to the simulation, blocking the simulation.
	 *
	 * This function blocks the simulation and it is slow, so it should be avoided.
	 *
	 * The actual number of particles pushed may be different than requested if
	 * the buffer is still (partially) full or smaller than requested.
	 *
	 * \param particles Pointer to an array of particles to be added.
	 * \param tags      Pointer to the corresponding array of tags.
	 * \param N         Number of particles to be added.
	 *
	 * \return Number of particles added to the simulation.
	 */
	inline CPU particle_index_t push(particle* particles, uint32_t* tags, particle_index_t N);

	/**
	 * \brief Set the detected particles to terminated, calling a callback
	 *        before doing so.
	 *
	 * This function blocks the simulation and it is slow, so it should be avoided.
	 *
	 * The callback function receives the particle data and the tag that
	 * belonged to the primary electron that initated the cascade the detected
	 * particle is part of.
	 *
	 * \param function Callback function to be called for each detected particle.
	 *                 Should have signature void(particle const &, uint32_t).
	 */
	template<typename detect_function>
	CPU void flush_detected(detect_function function);

	/**
	 * \brief Get number of particles currently in the simulation.
	 *
	 * More specifically, those that are not terminated or detected.
	 *
	 * This function blocks the simulation and it is slow, so it should be avoided.
	 */
	CPU particle_index_t get_running_count() const;

	/**
	 * \brief Get number of detected particles currently in the simulation.
	 *
	 * This function blocks the simulation and it is slow, so it should be avoided.
	 */
	CPU particle_index_t get_detected_count() const;



	/**
	 * \brief Allocate "input buffers" for asynchronously pushing new particles
	 *        to the simulation.
	 *
	 * \param N Size of the buffer (maximum number particles pushed).
	 *
	 * \see push_to_buffer
	 * \see push_to_simulation
	 */
	CPU void allocate_input_buffers(particle_index_t N);

	/**
	 * \brief Push particles to the input buffer, not adding them to the simulation.
	 *
	 * This function may be run asynchronously to the simulation.
	 * The actual number of particles pushed may be different than requested if
	 * the buffer is still (partially) full or smaller than requested.
	 *
	 * \param pool The work pool to add particles from.
	 *
	 * \see allocate_input_buffers
	 * \see push_to_simulation
	 */
	CPU void push_to_buffer(work_pool& pool);

	/**
	 * \brief Push electrons in the input buffer to the simulation.
	 *
	 * \see allocate_input_buffers
	 * \see push_to_buffer
	 */
	CPU void push_to_simulation();


	/**
	 * \brief Buffer the detected electrons in the current simulation memory to
	 *        an "output buffer", and free their memory for the main simulation.
	 *
	 * Technically, this function copies the entire simulation state into a
	 * buffer on the GPU (very fast). Other functions then copy this buffer to
	 * the CPU.
	 *
	 * \see flush_buffered
	 */
	CPU void buffer_detected();

	/**
	 * \brief Call callback function for each detected electron in the output
	 *        buffer.
	 *
	 * This function may be run asynchronously to the simulation.
	 *
	 * The callback function receives the particle data and the tag that
	 * belonged to the primary electron that initated the cascade the detected
	 * particle is part of.
	 *
	 * \param function Callback function to be called for each detected particle.
	 *                 Should have signature void(particle const &, uint32_t).
	 * \return Number of active electrons still in the simulation.
	 *
	 * \see buffer_detected
	 */
	template<typename detect_function>
	CPU particle_index_t flush_buffered(detect_function function);

	/// Perform a single iteration of the simulation for all particles.
	CPU void do_iteration();

	// Energy threshold (w.r.t. the vacuum level) below which particles must be terminated.
	real energy_threshold;

private:
	particle_manager_t _particles;
	material_manager_t _materials;
	geometry_manager_t _geometry;
	intersect_t _intersect;

	// GPU run settings
	unsigned int _threads_per_block = 32;
	unsigned int _num_blocks = 0;

	// Random number generators
	util::random_generator<true>* curand_states = nullptr;

	// Functions called by do_iteration()
	inline CPU void init();
	inline CPU void events();


	// Stream used for copying to/from input/output buffers.
	cudaStream_t buffer_stream = {};

	// Input buffers. "din" are on device, "hin" on host.
	particle_index_t buffer_in_size       = 0;
	bool*            buffer_din_data      = nullptr; // Value 0: no new electron, value 1: new electron
	particle*        buffer_din_particles = nullptr;
	uint32_t*        buffer_din_tags      = nullptr;
	bool*            buffer_hin_data      = nullptr;
	particle*        buffer_hin_particles = nullptr;
	uint32_t*        buffer_hin_tags      = nullptr;

	// Output buffers. "dout" are on device, "hout" on host.
	status_t* buffer_dout_status    = nullptr;
	particle* buffer_dout_particles = nullptr;
	uint32_t* buffer_dout_tags      = nullptr;
	status_t* buffer_hout_status    = nullptr;
	particle* buffer_hout_particles = nullptr;
	uint32_t* buffer_hout_tags      = nullptr;

	gpu_driver(gpu_driver const &) = delete;
	gpu_driver& operator=(gpu_driver const &) = delete;
};

}} // namespace nbl::drivers

#include "gpu_driver.inl"

#endif // __GPU_DRIVER_H_
