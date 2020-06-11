namespace nbl { namespace drivers {

template<typename material_t, bool gpu_flag>
CPU material_manager<material_t, gpu_flag>
	material_manager<material_t, gpu_flag>::create(std::vector<material_t> const & material_vector)
{
	return detail::material_manager_factory<material_t, gpu_flag>::create(material_vector);
}

template<typename material_t, bool gpu_flag>
CPU void material_manager<material_t, gpu_flag>::destroy(material_manager<material_t, gpu_flag> & manager)
{
	detail::material_manager_factory<material_t, gpu_flag>::free(manager);
}

template<typename material_t, bool gpu_flag>
PHYSICS material_t & material_manager<material_t, gpu_flag>::operator[](material_index_t i)
{
	return _materials[i];
}
template<typename material_t, bool gpu_flag>
PHYSICS material_t const & material_manager<material_t, gpu_flag>::operator[](material_index_t i) const
{
	return _materials[i];
}

template<typename material_t, bool gpu_flag>
PHYSICS bool material_manager<material_t, gpu_flag>::is_physical(material_index_t material_idx) const
{
	return material_idx >= 0;
}

namespace detail
{
	template<typename material_t>
	struct material_manager_factory<material_t, false>
	{
		using material_manager_t = material_manager<material_t, false>;
		using material_index_t = typename material_manager_t::material_index_t;

		inline static CPU material_manager_t create(std::vector<material_t> const & material_vector)
		{
			material_manager_t manager;

			if (material_vector.size() > std::numeric_limits<material_index_t>::max())
				throw std::runtime_error("Too many materials submitted");
			manager._capacity = static_cast<material_index_t>(material_vector.size());

			manager._materials = reinterpret_cast<material_t*>(malloc(manager._capacity*sizeof(material_t)));
			for (material_index_t i = 0; i < manager._capacity; ++i)
				manager._materials[i] = material_vector[i];

			return manager;
		}

		inline static CPU void free(material_manager_t & manager)
		{
			::free(manager._materials);
			manager._materials = nullptr;
			manager._capacity = 0;
		}
	};

#if CUDA_COMPILER_AVAILABLE
	template<typename material_t>
	struct material_manager_factory<material_t, true>
	{
		using material_manager_t = material_manager<material_t, true>;
		using material_index_t = typename material_manager_t::material_index_t;

		inline static CPU material_manager_t create(std::vector<material_t> const & material_vector)
		{
			material_manager_t manager;

			if (material_vector.size() > std::numeric_limits<material_index_t>::max())
				throw std::runtime_error("Too many materials submitted");
			manager._capacity = static_cast<material_index_t>(material_vector.size());

			cuda::cuda_new<material_t>(&manager._materials, manager._capacity);
			cuda::cuda_mem_scope<material_t>(manager._materials, manager._capacity,
				[&material_vector](material_t* device)
				{
					for (size_t i = 0; i < material_vector.size(); ++i)
						device[i] = material_vector[i];
				});

			return manager;
		}

		inline static CPU void free(material_manager_t & manager)
		{
			cudaFree(manager._materials);
			manager._materials = nullptr;
			manager._capacity = 0;
		}
	};
#endif // CUDA_COMPILER_AVAILABLE
} // namespace detail

}} // namespace nbl::drivers
