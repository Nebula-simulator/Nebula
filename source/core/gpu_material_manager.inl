namespace nbl {

template<typename material_t>
template<typename source_material_t>
CPU gpu_material_manager<material_t> gpu_material_manager<material_t>::create(
	cpu_material_manager<source_material_t> const & source)
{
	gpu_material_manager target;

	// Convert materials to GPU-compatible versions
	std::vector<material_t> material_vector;
	for (material_index_t i = 0; i < source.size(); ++i)
		material_vector.push_back(material_t::create(source[i]));

	// Copy all other data
	target._capacity = source.size();
	cuda::cuda_new<material_t>(&target._materials, target._capacity);
	cuda::cuda_mem_scope<material_t>(target._materials, target._capacity,
		[&material_vector](material_t* device)
		{
			for (size_t i = 0; i < material_vector.size(); ++i)
				device[i] = material_vector[i];
		});

	return target;
}

template<typename material_t>
CPU void gpu_material_manager<material_t>::destroy(gpu_material_manager<material_t> & manager)
{
	// Destroy all material data
	cuda::cuda_mem_scope<material_t>(manager._materials, manager._capacity,
		[&manager](material_t* materials){
			for (material_index_t i = 0; i < manager._capacity; ++i)
				material_t::destroy(materials[i]);
		});

	// And the rest
	cudaFree(manager._materials);
	manager._materials = nullptr;
	manager._capacity = 0;
}

template<typename material_t>
PHYSICS material_t & gpu_material_manager<material_t>::operator[](material_index_t i)
{
	return _materials[i];
}
template<typename material_t>
PHYSICS material_t const & gpu_material_manager<material_t>::operator[](material_index_t i) const
{
	return _materials[i];
}

template<typename material_t>
PHYSICS bool gpu_material_manager<material_t>::is_physical(material_index_t material_idx) const
{
	return material_idx >= 0;
}

} // namespace nbl
