#ifndef __GPU_PARTICLE_MANAGER_H_
#define __GPU_PARTICLE_MANAGER_H_

namespace nbl { namespace drivers {

template<typename material_manager_t>
class gpu_particle_manager
{
public:
	using status_t = uint8_t;
	using particle_index_t = uint32_t;
	using material_index_t = typename material_manager_t::material_index_t;

	static CPU gpu_particle_manager create(particle_index_t capacity);
	static CPU void destroy(gpu_particle_manager & manager);

	// Push particles to GPU. Returns how many particles were actually pushed.
	CPU particle_index_t push(particle* particles, uint32_t* tags, particle_index_t N);

	template<typename detect_function>
	CPU void flush_detected(detect_function function);

	CPU particle_index_t get_running_count() const;
	CPU particle_index_t get_detected_count() const;

	// Sort particles by event type. Used to prevent warp divergence.
	CPU void sort();

	// Direct access, no bounds checking
	inline PHYSICS particle & operator[](particle_index_t i);
	inline PHYSICS particle const & operator[](particle_index_t i) const;

	// Get max number of particles we can hold
	inline PHYSICS particle_index_t get_capacity() const;
	// Is there a memory location for this particle?
	inline PHYSICS bool exists(particle_index_t i) const;
	// Is particle active == not PENDING, DETECTED or TERMINATED
	inline PHYSICS bool active(particle_index_t i) const;

	// Get sorted particle index
	inline PHYSICS particle_index_t get_particle_index(particle_index_t i) const;
	// Get current material
	inline PHYSICS material_index_t get_material_index(particle_index_t i) const;
	inline PHYSICS void set_material_index(particle_index_t particle_idx, material_index_t new_material_idx);
	// Get last intersected triangle for a particle (or nullptr)
	inline PHYSICS triangle const * get_last_triangle(particle_index_t i) const;
	inline PHYSICS void forget_last_triangle(particle_index_t i);

	// Is next event elastic / inelastic / intersect?
	inline PHYSICS bool next_elastic(particle_index_t i) const;
	inline PHYSICS bool next_inelastic(particle_index_t i) const;
	inline PHYSICS bool next_intersect(particle_index_t i) const;
	inline PHYSICS bool is_detected(particle_index_t i) const;
	inline PHYSICS bool is_terminated(particle_index_t i) const;

	// Add a new particle, without checking that the target index is actually free.
	inline PHYSICS void add_particle(particle_index_t target_idx, particle new_particle, uint32_t new_tag);
	// Is a secondary slot free for this thread to write in?
	inline PHYSICS bool secondary_slot_free() const;
	// Assumes that secondary_slot_free() returns true!
	inline PHYSICS void create_secondary(particle_index_t primary_idx, particle secondary_particle);

	// Terminate particle
	inline PHYSICS void terminate(particle_index_t i);
	// Detect particle
	inline PHYSICS void detect(particle_index_t i);
	// Set next scattering event
	inline PHYSICS void set_scatter_event(particle_index_t i, uint8_t event);
	inline PHYSICS void set_intersect_event(particle_index_t i, triangle* t);

	// Set particle to pending (for an inelastic event)
	inline PHYSICS void pending(particle_index_t i);
	// Set particle to inelastic (for previously pending event)
	inline PHYSICS void inelastic(particle_index_t i);

private:
	inline PHYSICS particle_index_t get_secondary_slot() const;

	particle_index_t  _capacity      = 0;       // How many particles can we hold
	status_t*         _status        = nullptr; // status_enum for each particle
	particle_index_t* _particle_idx  = nullptr; // Particle index, sorted by status
	particle*         _particles     = nullptr; // The actual particle data (position, direction, energy)
	uint32_t*         _tags          = nullptr; // Each particle has an associated tag
	material_index_t* _material_idx  = nullptr; // Current material the particle is in
	triangle**        _last_triangle = nullptr; // Pointer to last intersected triangle

	// Sorting stuff
	particle_index_t* _radix_index = nullptr;
	status_t* _radix_dump          = nullptr;
	void* _radix_temp              = nullptr;
	size_t _radix_temp_size        = 0;

	enum status_enum : status_t
	{
		PENDING         = 0b0000,
		INELASTIC_EVENT = 0b0100,
		ELASTIC_EVENT   = 0b0001,
		INTERSECT_EVENT = 0b0010,
		NO_EVENT        = 0b0110,
		DETECTED        = 0b1010,
		NEW_SECONDARY   = 0b1110,
		TERMINATED      = 0b0011
	};

	template<typename scatter_list_t,
		typename intersect_t,
		typename geometry_manager_t
	>
	friend class gpu_driver;
};

}} // namespace nbl::drivers

#include "gpu_particle_manager.inl"

#endif // __GPU_PARTICLE_MANAGER_H_
