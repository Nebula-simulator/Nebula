#ifndef __MATERIAL_MANAGER_H_
#define __MATERIAL_MANAGER_H_

namespace nbl { namespace drivers {

namespace detail
{
	/**
	 * \brief Responsible for managing the memory on the CPU or GPU device.
	 */
	template<typename material_type, bool gpu_flag>
	struct material_manager_factory;
}

/**
 * \brief Provides access to all materials in the simulation.
 *
 * The material manager is allocated and deallocated by the static
 * {@link create} and {@link destroy} functions, there is no constructor or
 * destructor to be used.
 *
 * \tparam material_type Datatype of the material. Should be some specialization
 *                       of ::material.
 * \tparam gpu_flag      Set to true if simulation is to be run on the GPU,
 *                       false for CPU.
 */
template<typename material_type, bool gpu_flag>
class material_manager
{
public:
	using material_t = material_type;
	using material_index_t = int;     ///< Type used to index materials.

	/**
	 * \brief Allocate memory for the data, provided a list of materials.
	 */
	static CPU material_manager<material_t, gpu_flag> create(
		std::vector<material_t> const & materials);
	/**
	 * \brief Destroy the material manager, deallocating the data.
	 */
	static CPU void destroy(material_manager & manager);

	/**
	 * \brief Direct read-write access to a material.
	 *
	 * No bounds checking is done.
	 *
	 * \param i Index to the material.
	 */
	inline PHYSICS material_t & operator[](material_index_t i);
	/**
	 * \brief Direct read-only access to a material.
	 *
	 * No bounds checking is done.
	 *
	 * \param i Index to the material.
	 */
	inline PHYSICS material_t const & operator[](material_index_t i) const;

	/**
	 * \brief Checks if the material is not a "special" vacuum-type material.
	 *
	 * This function only checks if the provided index >= 0, that is, it does
	 * not check if the material designated by the index actually exists.
	 *
	 * \see index_enum
	 */
	inline PHYSICS bool is_physical(material_index_t material_idx) const;

	/**
	 * \brief Enumeration of "special" materials, such as vacuum and detectors.
	 *
	 * When an electron goes crosses a boundary from a material to one of these,
	 * something special happens according to the particular value in this enum.
	 */
	enum index_enum : material_index_t
	{
		NOP = -128,            ///< No-op. TODO: support this
		TERMINATOR = -127,     ///< Destroy the electron
		DETECTOR = -126,       ///< Flag as detected
		DETECTOR_LT50 = -125,  ///< Flag as detected, if energy < 50 eV
		DETECTOR_GE50 = -124,  ///< Flag as detected, if energy >= 50 eV
		VACUUM = -123,         ///< Go to vacuum propagation
		MIRROR = -122          ///< Specularly reflect
	};

private:
	material_index_t _capacity = 0;
	material_t* _materials = nullptr;

	friend struct detail::material_manager_factory<material_t, gpu_flag>;
};

}} // namespace nbl::drivers

#include "material_manager.inl"

#endif // __MATERIAL_MANAGER_H_
