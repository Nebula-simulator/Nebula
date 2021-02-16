#ifndef __GPU_MATERIAL_MANAGER_H_
#define __GPU_MATERIAL_MANAGER_H_

#include "cpu_material_manager.h"

namespace nbl {

/**
 * \brief Provides access to all materials in the simulation.
 *
 * The material manager is allocated and deallocated by the static
 * {@link create} and {@link destroy} functions, there is no constructor or
 * destructor. Reason is that the copy constructor is called whenever you pass
 * a class to a CUDA kernel.
 *
 * \tparam material_type Datatype of the material. Should be some specialization
 *                       of ::material.
 */
template<typename material_type>
class gpu_material_manager
{
public:
	using material_t = material_type;
	using material_index_t = int;     ///< Type used to index materials.

	/**
	 * \brief Allocate memory for the data.
	 *
	 * This function clones data from a CPU material manager. This "source"
	 * material manager may have a different material type than the GPU material
	 * manager. This is to allow conversion from "CPU" materials to "GPU"
	 * materials. Usually, this conversion means that large data sets are copied
	 * to the GPU.
	 *
	 * \param materials Material manager to copy from
	 */
	template<typename source_material_t>
	static CPU gpu_material_manager<material_t> create(
		cpu_material_manager<source_material_t> const & materials);

	/**
	 * \brief Destroy the material manager, deallocating the data.
	 */
	static CPU void destroy(gpu_material_manager & manager);

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
	inline PHYSICS bool is_physical(material_index_t material_idx) const;

private:
	material_index_t _capacity = 0;
	material_t* _materials = nullptr;
};

} // namespace nbl

#include "gpu_material_manager.inl"

#endif // __GPU_MATERIAL_MANAGER_H_
