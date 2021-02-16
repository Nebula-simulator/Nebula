#include <initializer_list>

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

	template<typename... scatter_types, size_t... is>
	void destroy_scatter_list(
		scatter_list<scatter_types...> & list,
		nbl::index_sequence<is...>)
	{
		(void)std::initializer_list<int>{
			(scatter_types::destroy(nbl::tuple::get<is>(list)), 0)... };
	}
}

template<typename... scatter_types>
CPU auto material<scatter_list<scatter_types...>>::create(nbl::hdf5_file const & mat)
	-> material
{
	material<base_t> target;

	target.barrier = static_cast<real>(mat.get_property_quantity("barrier") / nbl::units::eV);
	target.base_t::operator=(base_t(scatter_types::create(mat)...));

	return target;
}

template<typename... scatter_types>
template<typename... source_scatter_types>
CPU auto material<scatter_list<scatter_types...>>::create(
	material<scatter_list<source_scatter_types...>> const & source)
	-> material
{
	material<base_t> target;

	target.barrier = source.barrier;
	target.base_t::operator=(detail::create_scatter_list<
		scatter_list<source_scatter_types...>,
		scatter_types...>(source, nbl::make_index_sequence<sizeof...(scatter_types)>{}));

	return target;
}

template<typename... scatter_types>
CPU void material<scatter_list<scatter_types...>>::destroy(material & mat)
{
	detail::destroy_scatter_list(mat,
		nbl::make_index_sequence<sizeof...(scatter_types)>{});
}
