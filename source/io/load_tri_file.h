#ifndef __LOAD_TRI_FILE_H_
#define __LOAD_TRI_FILE_H_

#include <vector>
#include <string>
#include "../core/triangle.h"

namespace nbl {

	/**
	 * \brief Load triangle file.
	 *
	 * This is a text file, where each triangle is stored as a single line.
	 * The first two numbers on a line are integers representing the material
	 * indices on either side of the triangle. The other 9 numbers represent the
	 * x,y,z coordinates (in nm) of the triangle vertices.
	 *
	 * \return A vector of the triangles.
	 */
	std::vector<triangle> load_tri_file(std::string const & filename);

} // namespace nbl

#endif // __LOAD_TRI_FILE_H_
