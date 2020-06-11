namespace detail
{
	template<typename source_type, typename... target_scatter_types, size_t... is>
	scatter_list<target_scatter_types...> create_scatter_list(
		source_type const & source,
		nbl::index_sequence<is...>)
	{
		return { target_scatter_types::create(
			nbl::tuple::get<is>(source))... };
	}
}

template<typename... scatter_types>
CPU material<scatter_list<scatter_types...>>::material(nbl::hdf5_file const & mat)
	: scatter_list<scatter_types...>(scatter_types::create(mat)...),
	barrier(static_cast<real>(mat.get_property_quantity("barrier") / nbl::units::eV))
{
}

template<typename... scatter_types>
template<typename... source_scatter_types>
CPU material<scatter_list<scatter_types...>>::material(
	material<scatter_list<source_scatter_types...>> const & source)
	:
	scatter_list<scatter_types...>(detail::create_scatter_list<
		scatter_list<source_scatter_types...>,
		scatter_types...>(source, nbl::make_index_sequence<sizeof...(scatter_types)>{})),
	barrier(source.barrier)
{
}
