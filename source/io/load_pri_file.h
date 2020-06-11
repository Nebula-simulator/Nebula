#ifndef __LOAD_PRI_FILE_H_
#define __LOAD_PRI_FILE_H_

#include <vector>
#include <tuple>
#include <string>
#include "../core/particle.h"

#if !CUDA_HEADERS_AVAILABLE
struct int2 { int x, y; };
#endif

namespace nbl {

	/**
	 * \brief Load primary electrons file.
	 *
	 * This is a binary file, where each electron is stored as 7 floats and 2 ints.
	 * The values are
	 *  - The position x,y,z in nm
	 *  - Direction vector x,y,z (unnormalised)
	 *  - Energy in eV
	 *  - Pixel indices x and y
	 *
	 * \return A vector of the particle data,
	 *         and a corresponding vector of pixel indices.
	 */
	std::pair<std::vector<particle>, std::vector<int2>> load_pri_file(
		std::string const & filename,
		vec3 min_pos,
		vec3 max_pos, real max_E);

	/**
	 * \brief Sort primary electrons by their starting positions.
	 *
	 * Sorting is done by Morton index, to make sure that nearby particles in
	 * the exposure are also spatially close. This gives a significant speedup
	 * in the GPU collision detection routine, but in practice, the gains are
	 * usually less than the cost of sorting, especially if the input is already
	 * partially sorted.
	 */
	void sort_pri_file(
		std::vector<particle>& particle_vec,
		std::vector<int2>& pixel_vec);

	/**
	 * \brief Shuffle the primary particles such that the first `prescan_size`
	 *        are uniformly sampled.
	 *
	 * The others are left mostly untouched.
	 */
	void prescan_shuffle(
		std::vector<particle>& particle_vec,
		std::vector<int2>& pixel_vec,
		size_t prescan_size);

} // namespace nbl

#endif // __LOAD_PRI_FILE_H_
