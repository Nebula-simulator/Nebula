#ifndef __CPU_MATERIAL_MANAGER_H_
#define __CPU_MATERIAL_MANAGER_H_

namespace nbl {

/**
 * \brief Provides access to all materials in the simulation.
 *
 * \tparam material_type Datatype of the material. Should be some specialization
 *                       of ::material.
 */
template<typename material_type>
class cpu_material_manager
{
public:
	using material_t = material_type;
	using material_index_t = int;     ///< Type used to index materials.

	cpu_material_manager() = default;
	~cpu_material_manager();

	/**
	 * \brief Load a new material from a HDF5 file and add it to the list.
	 *
	 * \param mat The HDF5 file to load the data from
	 */
	void add(hdf5_file const & mat);

	/**
	 * \brief Get the number of materials.
	 */
	material_index_t size() const;

	/**
	 * \brief Get the maximal electron energy accepted by all materials.
	 */
	real get_max_energy() const;

	/**
	 * \brief Direct read-write access to a material.
	 *
	 * No bounds checking is done.
	 *
	 * \param i Index to the material.
	 */
	inline material_t & operator[](material_index_t i);
	/**
	 * \brief Direct read-only access to a material.
	 *
	 * No bounds checking is done.
	 *
	 * \param i Index to the material.
	 */
	inline material_t const & operator[](material_index_t i) const;

	/**
	 * \brief Checks if the material is not a "special" vacuum or detector
	 *        material.
	 *
	 * This function only checks if the provided index >= 0, that is, it does
	 * not check if the material designated by the index actually exists.
	 *
	 * \param material_idx Material index to check.
	 *
	 * \see special_materials
	 */
	inline bool is_physical(material_index_t material_idx) const;

private:
	std::vector<material_t> _materials;

	cpu_material_manager(cpu_material_manager&) = delete;
	cpu_material_manager& operator=(cpu_material_manager const &) = delete;
};

/**
 * \brief Enumeration of "special" materials, such as vacuum and detectors.
 *
 * When an electron goes crosses a boundary from a material to one of these,
 * something special happens according to the particular value in this enum.
 */
enum special_materials
{
	NOP = -128,            ///< No-op.
	TERMINATOR = -127,     ///< Destroy the electron
	DETECTOR = -126,       ///< Flag as detected
	DETECTOR_LT50 = -125,  ///< Flag as detected, if energy < 50 eV
	DETECTOR_GE50 = -124,  ///< Flag as detected, if energy >= 50 eV
	VACUUM = -123,         ///< Go to vacuum propagation
	MIRROR = -122          ///< Specularly reflect
};

} // namespace nbl

#include "cpu_material_manager.inl"

#endif // __CPU_MATERIAL_MANAGER_H_
