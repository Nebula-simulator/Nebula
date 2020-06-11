#ifndef __WORK_POOL_H_
#define __WORK_POOL_H_

#include <mutex>
#include "../core/particle.h"

/**
 * \brief Very simple class to keep track of work to be done in a thread-safe way.
 *
 * This class does not take ownership of the data.
 * It simply assumes that the particles are not deleted until the work is done.
 */

class work_pool
{
public:
	/**
	 * \brief Construct to invalid state.
	 */
	explicit work_pool()
		: next_primary(nullptr), next_tag(nullptr), primaries_to_go(0)
	{}

	/**
	 * \brief Constructor.
	 *
	 * \param primaries Pointer to primary electron data
	 * \param tags      Pointer to tag data
	 * \param N         Number of particles to be simulated
	 */
	work_pool(particle* primaries, uint32_t* tags, size_t N) :
		next_primary(primaries), next_tag(tags), primaries_to_go(N)
	{}

	/**
	 * \brief Get work to be done
	 *
	 * \param batch_size Number of particles requested
	 *
	 * \return Tuple of
	 *   - Pointer to first particle
	 *   - Pointer to first corresponding tag
	 *   - Number of particles obtained, may be less than batch_size
	 */
	std::tuple<particle*, uint32_t*, uint32_t> get_work(uint32_t batch_size)
	{
		std::lock_guard<std::mutex> lock(mutex);

		uint32_t particles_pushed = (uint32_t)std::min<size_t>(batch_size, primaries_to_go);
		auto return_data = std::make_tuple(next_primary, next_tag, particles_pushed);
		next_primary += particles_pushed;
		next_tag += particles_pushed;
		primaries_to_go -= particles_pushed;

		return return_data;
	}

	/**
	 * Get work to be done
	 *
	 * This function checks `in_data` to see whether a slot is free. The
	 * particle and tag data are then copied, and `in_data` is updated.
	 *
	 * \param data      Indicates if a slot is free (0 if free, 1 if full).
	 * \param particles Array for particle data
	 * \param tags      Array for tag data
	 * \param N         Length of the arrays
	 */
	void get_work(bool* data, particle* particles, uint32_t* tags, uint32_t N)
	{
		std::lock_guard<std::mutex> lock(mutex);

		for (uint32_t i = 0; i < N; ++i)
		{
			if (primaries_to_go == 0)
				break;
			if (data[i] != 0)
				continue;

			data[i] = 1;
			particles[i] = *next_primary;
			tags[i] = *next_tag;

			++next_primary;
			++next_tag;
			--primaries_to_go;
		}
	}

	/// Get the amount of work left in the pool
	size_t get_primaries_to_go() const
	{
		std::lock_guard<std::mutex> lock(mutex);
		return primaries_to_go;
	}

	/// Test if the amount of work is equal to zero.
	bool done() const
	{
		return get_primaries_to_go() == 0;
	}

	work_pool& operator=(work_pool const & rhs)
	{
		if (this == &rhs)
			return *this;

		std::lock_guard<std::mutex> lock1(rhs.mutex);
		std::lock_guard<std::mutex> lock2(mutex);
		next_primary = rhs.next_primary;
		next_tag = rhs.next_tag;
		primaries_to_go = rhs.primaries_to_go;

		return *this;
	}

private:
	mutable std::mutex mutex;

	particle* next_primary;
	uint32_t* next_tag;
	size_t primaries_to_go;
};

#endif // __WORK_POOL_H_
