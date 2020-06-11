#ifndef __RANGE_H_
#define __RANGE_H_

#include <cstddef>
#include <cmath>
#include <utility>

namespace nbl { namespace util {

/**
 * \brief CRTP base class for ranges
 *
 * \tparam concrete_range Derived concrete range type
 */
template<class concrete_range>
class basic_range
{
private:
	class const_iterator_impl
	{
	public:
		using value_type = typename concrete_range::value_type;
		using size_type = typename concrete_range::size_type;
		using difference_type = typename concrete_range::difference_type;

		const_iterator_impl(concrete_range const & range, size_type index)
			: r(range), i(index)
		{}

		value_type operator*() const
		{
			return r[i];
		}

		bool operator==(const const_iterator_impl& rhs) const
		{
			return (i == rhs.i && r == rhs.r);
		}
		bool operator!=(const const_iterator_impl& rhs) const
		{
			return !(*this == rhs);
		}

		const_iterator_impl& operator++()
		{
			++i;
			return *this;
		}
		const_iterator_impl operator++(int)
		{
			const_iterator_impl temp = *this;
			++i;
			return temp;
		}

		const_iterator_impl& operator--()
		{
			--i;
			return *this;
		}
		const_iterator_impl operator--(int)
		{
			const_iterator_impl temp = *this;
			--i;
			return temp;
		}

	private:
		concrete_range r;
		size_type i;
	};

public:
	// Member types
	using size_type = std::size_t;              ///< Indexing type
	using difference_type = std::ptrdiff_t;     ///< Difference type
	using iterator = const_iterator_impl;       ///< Random access iterator
	using const_iterator = const_iterator_impl; ///< Constant iterator
	// using reverse_iterator =
	// using const_reverse_iterator =

	/// Get iterator to start of range
	const_iterator begin() const
	{
		return const_iterator(*static_cast<concrete_range const *>(this), 0);
	}
	/// Get iterator to one-past end of range
	const_iterator end() const
	{
		return const_iterator(*static_cast<concrete_range const *>(this), count);
	}

	/// Get number of elements in the range
	size_type size() const
	{
		return count;
	}
	/// Whether the range is empty (equivalent to size() == 0)
	bool empty() const
	{
		return count == 0;
	}

protected:
	explicit basic_range(size_type count) :
		count(count)
	{}

	/// Equality operator
	bool operator==(basic_range const & rhs) const
	{
		return count == rhs.count;
	}

	size_type count;
};


/**
 * \brief Range of evenly-spaced integers
 */
template<typename T = int>
class range : public basic_range<range<T>>
{
	using base_t = basic_range<range<T>>;

public:
	using value_type = T;
	using typename base_t::size_type;

	/**
	 * \brief Construct range from start to stop with a given step size.
	 *
	 * Note: \p stop is not included.
	 *
	 * \param start First value in the range
	 * \param stop  Stopping value in the range (not included)
	 * \param step  Step size
	 */
	range(value_type start, value_type stop, value_type step) :
		first(start),
		step(step),
		base_t((stop-start)/step)
	{}

	/**
	 * \brief Construct range from start to stop. Step size is +-1.
	 *
	 * The sign of the step size is determined automatically.
	 * Note: \p stop is not included.
	 *
	 * \param start First value in the range
	 * \param stop  Stopping value in the range (not included)
	 */
	range(value_type start, value_type stop) :
		first(start),
		step(stop>start ? 1 : -1),
		base_t(stop>start ? stop-start : start-stop)
	{}
	/**
	 * \brief Construct range from 0 to stop. Step size is +-1.
	 *
	 * The sign of the step size is determined automatically.
	 * Note: \p stop is not included.
	 *
	 * \param stop  Stopping value in the range (not included)
	 */
	range(value_type stop) :
		first(0),
		step(stop>0 ? 1 : -1),
		base_t(stop>0 ? stop : -stop)
	{}

	/// Read-only access first value
	value_type front() const
	{
		return first;
	}
	/// Read-only access last value
	value_type back() const
	{
		return (*this)[base_t::count-1];
	}

	/// Random access to value
	value_type operator[](size_type index) const
	{
		return first + index * step;
	}

	/// Equality comparison
	bool operator==(range const & rhs) const
	{
		return base_t::operator==(rhs)
			&& first == rhs.first
			&& step == rhs.step;
	}

private:
	value_type first;
	value_type step;
};

/**
 * \brief Range of evenly-spaced real numbers.
 *
 * Due to numerical round-off errors, it is not guaranteed that the last value
 * in the range is exactly equal to what was provided.
 */
template<typename T = double>
class linspace : public basic_range<linspace<T>>
{
	using base_t = basic_range<linspace<T>>;

public:
	using value_type = T;
	using typename base_t::size_type;

	/**
	 * \brief Construct range from start to stop in a given number of steps.
	 *
	 * \p start and \p stop are both included in the range.
	 *
	 * \param start First value in the range
	 * \param stop  Last value in the range
	 * \param num   Number of values in the range
	 */
	linspace(value_type start, value_type stop, size_type num) :
		base_t(num),
		first(start),
		step((stop-start)/(num-1))
	{}

	/// Read-only access first value
	value_type front() const
	{
		return first;
	}
	/// Read-only access last value
	value_type back() const
	{
		return (*this)[base_t::count-1];
	}

	/// Random access to value
	value_type operator[](size_type index) const
	{
		return first + index * step;
	}

	/// Equality comparison
	bool operator==(linspace const & rhs) const
	{
		return base_t::operator==(rhs)
			&& first == rhs.first
			&& step == rhs.step;
	}

private:
	value_type first;
	value_type step;
};

/**
 * \brief Range of logarithmically-spaced real numbers.
 *
 * Due to numerical round-off errors, it is not guaranteed that the last value
 * in the range is exactly equal to what was provided.
 */
template<typename T = double>
class geomspace : public basic_range<geomspace<T>>
{
	using base_t = basic_range<geomspace<T>>;

public:
	using value_type = T;
	using typename base_t::size_type;

	/**
	 * \brief Construct range from start to stop in a given number of steps.
	 *
	 * \p start and \p stop are both included in the range.
	 *
	 * \param start First value in the range
	 * \param stop  Last value in the range
	 * \param num   Number of values in the range
	 */
	geomspace(value_type start, value_type stop, size_type num) :
		base_t(num),
		first(start),
		step(std::log(stop/start) / (num - 1))
	{}

	/// Read-only access first value
	value_type front() const
	{
		return first;
	}
	/// Read-only access last value
	value_type back() const
	{
		return (*this)[base_t::count-1];
	}

	/// Random access to value
	value_type operator[](size_type index) const
	{
		return first * std::exp(index * step);
	}

	/// Equality comparison
	bool operator==(geomspace const & rhs) const
	{
		return base_t::operator==(rhs)
			&& first == rhs.first
			&& step == rhs.step;
	}

private:
	value_type first;
	using step_type = decltype(std::declval<T>() / std::declval<T>());
	step_type step;
};

}} // namespace nbl::util

#endif // __RANGE_H_
