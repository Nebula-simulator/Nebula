#ifndef __MATERIAL_H_
#define __MATERIAL_H_

#include "scatter_list.h"
#include "../io/hdf5_file.h"

template<typename...>
struct material;

/**
 * \brief Material class.
 *
 * Template parameter is a ::scatter_list, i.e. a list of physical scattering
 * mechanisms. These scattering mechansims contain the physics of what actually
 * goes on inside the material.
 */
template<typename... scatter_types>
struct material<scatter_list<scatter_types...>>
	: public scatter_list<scatter_types...>
{
public:
	using base_t = scatter_list<scatter_types...>;

	/**
	 * \brief Read from a HDF5 material file.
	 *
	 * Under the hood, this calls the create() function for each of the
	 * \p scatter_types that are part of this material.
	 *
	 * \param mat HDF5 file with the data
	 */
	static CPU material create(nbl::hdf5_file const & mat);

	/**
	 * \brief Clone from a different material.
	 *
	 * The "source" material may have a different list of scattering types, as
	 * long as they are convertible. This happens, for example, when cloning
	 * a CPU material to a GPU.
	 *
	 * \param source Material class with data to clone
	 */
	template<typename... source_scatter_types>
	static CPU material create(material<scatter_list<source_scatter_types...>> const & source);

	/**
	 * \brief Deallocate all data
	 */
	static CPU void destroy(material & mat);

	/**
	 * \brief Barrier energy, that is, the work function plus the Fermi energy.
	 *
	 * An electron cannot escape the material unless it has more kinetic energy
	 * than this value.
	 */
	real barrier = 0;
};

#include "material.inl"

#endif // __MATERIAL_H_
