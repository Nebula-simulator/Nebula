#ifndef __TUPLE_H_
#define __TUPLE_H_

/**
 * \namespace nbl::tuple
 * \brief Contains something vaguely similar to `std::tuple`.
 *
 * We need something that works on both CUDA devices and host, for holding a
 * list of scattering mechanisms.
 *
 * NOTE THAT THIS IS NOT BY ANY MEANS A FULLY COMPLIANT TUPLE IMPLEMENTATION.
 * It cannot hold references, move-only types, or any of that.
 *
 * There is no perfect forwarding anywhere.
 */

#include "variadic.h"

namespace nbl { namespace tuple {

namespace detail
{
	template<size_t, typename T>
	struct tuple_element
	{
		tuple_element() = default;
		tuple_element(T value) : value_(value) {}
		T value_;
	};

	template<typename sequence, typename... types>
	struct tuple_impl;
	template<size_t... indices, typename... types>
	struct tuple_impl<index_sequence<indices...>, types...>
		: tuple_element<indices, types>...
	{
		explicit tuple_impl() = default;
		explicit tuple_impl(types... elements)
			: tuple_element<indices, types>(elements)...
		{}
	};

	#ifdef __NVCC__
	// TODO: Workaround for nvcc before CUDA 9
	template<size_t i, typename... types>
	struct ___get_helper { using type = type_at_index<i, types...>; };
	#endif
}


/**
 * \brief Tuple class, vaguely similar to `std::tuple`.
 *
 * This is a very basic implementation. It cannot hold things like references
 * or move-only types.
 */
template<typename... types>
struct tuple
///\cond INTERNAL # Hide the inheritance relationship in doxygen
	: detail::tuple_impl<make_index_sequence<sizeof...(types)>, types...>
///\endcond
{
	/**
	 * \brief The type this class inherits from
	 */
	using base_t = detail::tuple_impl<make_index_sequence<sizeof...(types)>, types...>;

	/**
	 * \brief Default constructor
	 */
	tuple() = default;

	/**
	 * \brief Construct from elements.
	 *
	 * No perfect forwarding.
	 */
	explicit tuple(types... elements) :
		base_t(elements...)
	{}
};


/**
 * \brief Get a reference to the `i`'th element in the tuple.
 *
 * \tparam i   Index to the element in the tuple
 * \param  tup Tuple to get element from
 */
template<size_t i, typename... types>
#ifdef __NVCC__
PHYSICS typename detail::___get_helper<i, types...>::type& get(tuple<types...>& tup)
#else
PHYSICS type_at_index<i, types...>& get(tuple<types...>& tup)
#endif
{
	detail::tuple_element<i, type_at_index<i, types...>>& base = tup;
	return base.value_;
}

/**
 * \brief Get a const reference to the `i`'th element in the tuple.
 *
 * \tparam i   Index to the element in the tuple
 * \param  tup Tuple to get element from
 */
template<size_t i, typename... types>
#ifdef __NVCC__
PHYSICS typename detail::___get_helper<i, types...>::type const & get(tuple<types...> const & tup)
#else
PHYSICS type_at_index<i, types...> const & get(tuple<types...> const & tup)
#endif
{
	detail::tuple_element<i, type_at_index<i, types...>> const & base = tup;
	return base.value_;
}


namespace detail
{
	// for_each()
	template<size_t i, typename T, typename F>
	inline PHYSICS typename std::enable_if<i != 0, void>::type for_each_impl(T const & tup, F& fun)
	{
		fun(get<i - 1>(tup), i - 1);
		for_each_impl<i - 1>(tup, fun);
	}
	template<size_t i, typename T, typename F>
	inline PHYSICS typename std::enable_if<i == 0, void>::type for_each_impl(T const &, F&)
	{}


	// visit_at()
	template<size_t i, typename T, typename F>
	inline PHYSICS typename std::enable_if<i != 0, void>::type visit_at_impl(T const & tup, size_t idx, F& fun)
	{
		if (idx == i - 1) fun(get<i - 1>(tup));
		else visit_at_impl<i - 1>(tup, idx, fun);
	}
	template<size_t i, typename T, typename F>
	inline PHYSICS typename std::enable_if<i == 0, void>::type visit_at_impl(T const &, size_t, F&)
	{}
}

/**
 * \brief Call a function on every element of a tuple.
 *
 * Currently, only a const version exists.
 *
 * \param tup The tuple.
 * \param fun The function to be called.
 */
template<typename F, typename... Ts>
PHYSICS void for_each(tuple<Ts...> const & tup, F& fun)
{
	detail::for_each_impl<sizeof...(Ts)>(tup, fun);
}

/**
 * \brief Call a function with a given element of the tuple.
 *
 * Similar to `fun(get<idx>(tup))`, with the exception that `idx` need not be
 * known at compile-time. If the index is known at compile-time, the other
 * method is recommended.
 *
 * Currently, only a const version exists.
 *
 * \param tup The tuple.
 * \param idx The index to the element of the tuple that the function needs to
 *            be called for.
 * \param fun The function to be called.
 */
template<typename F, typename... Ts>
PHYSICS void visit_at(tuple<Ts...> const & tup, size_t idx, F& fun)
{
	detail::visit_at_impl<sizeof...(Ts)>(tup, idx, fun);
}

}} // namespace nbl::tuple

#endif // __TUPLE_H_
