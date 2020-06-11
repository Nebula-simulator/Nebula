#ifndef __CPU_PARTICLE_MANAGER_H_
#define __CPU_PARTICLE_MANAGER_H_

#include <vector>
#include <map>

namespace nbl { namespace drivers {

/*
 * Generic CPU particle manager.
 * Serves as a base class for a "simple" CPU particle manager
 * and a version with extended tracking facilities.
 */
template<typename material_manager_t>
class cpu_particle_manager
{
public:
	using particle_index_t = size_t;
	using material_index_t = typename material_manager_t::material_index_t;
	using primary_tag_t = uint32_t;

	static cpu_particle_manager create();
	static void destroy(cpu_particle_manager & manager);

	// Push particles to GPU. Returns how many particles were actually pushed.
	particle_index_t push(particle* particles, primary_tag_t* tags, particle_index_t N);

	// Sets detected particles to terminated and calls callback.
	// Does not remove detected particles from memory.
	template<typename detect_function>
	void flush_detected(detect_function func);
	// Remove terminated particles from memory.
	void flush_terminated();

	particle_index_t get_total_count() const;
	particle_index_t get_running_count() const;
	particle_index_t get_detected_count() const;

	// Direct access, no bounds checking
	inline PHYSICS particle & operator[](particle_index_t i);
	inline PHYSICS particle const & operator[](particle_index_t i) const;

	// Is there a memory location for this particle?
	inline PHYSICS bool exists(particle_index_t i) const;
	// Is particle active == not DETECTED or TERMINATED
	inline PHYSICS bool active(particle_index_t i) const;

	// Get current material
	inline PHYSICS material_index_t get_material_index(particle_index_t i) const;
	inline PHYSICS void set_material_index(particle_index_t particle_idx, material_index_t new_material_idx);

	// Get primary tag
	inline PHYSICS primary_tag_t get_primary_tag(particle_index_t i) const;

	// Get last intersected triangle for a particle (or nullptr)
	inline PHYSICS triangle const * get_last_triangle(particle_index_t i) const;
	inline PHYSICS void forget_last_triangle(particle_index_t i);

	// Is next event scatter / intersect?
	inline PHYSICS bool next_scatter(particle_index_t i) const;
	inline PHYSICS uint8_t get_next_scatter(particle_index_t i) const;
	inline PHYSICS bool next_intersect(particle_index_t i) const;

	// Add new particle
	inline PHYSICS void create_secondary(particle_index_t primary_idx, particle secondary_particle);

	// Terminate particle
	inline PHYSICS void terminate(particle_index_t i);
	// Detect particle
	inline PHYSICS void detect(particle_index_t i);
	// Set next scattering event
	inline PHYSICS void set_scatter_event(particle_index_t i, scatter_event event);
	inline PHYSICS void set_intersect_event(particle_index_t i, intersect_event event);

protected:
	enum particle_status
	{
		SCATTER_EVENT,
		INTERSECT_EVENT,
		NO_EVENT,
		DETECTED,
		TERMINATED
	};

	/*
	 * particle_struct holds relevant data for each particle.
	 * It inherits from the "additional_data" template parameter, which may be
	 * void. This mechanism allows derived particle managers to store additional
	 * data in addition to the necessities.
	 */
	struct particle_struct
	{
		particle_status status;
		uint8_t next_scatter;
		material_index_t current_material;
		particle particle_data;
		primary_tag_t primary_tag; // Tag belonging to primary electron
		uint32_t secondary_tag;    // Unique tag for this electron in the primary's cascade
		triangle* last_triangle;
	};
	std::vector<particle_struct> particles;

	struct cascade_struct
	{
		uint32_t next_secondary_tag;
		uint32_t running_count;
	};
	std::map<primary_tag_t, cascade_struct> cascades;
};

}} // namespace nbl::drivers

#include "cpu_particle_manager.inl"

#endif // __CPU_PARTICLE_MANAGER_H_
