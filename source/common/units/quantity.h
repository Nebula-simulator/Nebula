#ifndef __QUANTITY_H_
#define __QUANTITY_H_

/*
 * A quantity is a number plus a dimension.
 */

#include <cmath>
#include <stdexcept>
#include "dimension.h"

namespace nbl { namespace units {

template<typename T>
struct quantity
{
	using value_type = T;

	value_type value;
	dimension units;

	explicit constexpr quantity(value_type val = value_type{},
		dimension unit = dimensions::dimensionless)
		: value(val), units(unit)
	{}

	bool dimensionless() const
	{
		return units == dimensions::dimensionless;
	}

	operator T() const
	{
		if (!dimensionless())
			throw std::runtime_error("Can only convert dimensionless quantity"
				" to scalar");
		return value;
	}

	quantity& operator+=(quantity const & rhs)
	{
		if (units != rhs.units)
			throw std::runtime_error("Adding incompatible units.");
		value += rhs.value;
		return *this;
	}
	quantity& operator-=(quantity const & rhs)
	{
		if (units != rhs.units)
			throw std::runtime_error("Adding incompatible units.");
		value -= rhs.value;
		return *this;
	}
	quantity& operator*=(quantity const & rhs)
	{
		value *= rhs.value;
		units *= rhs.units;
		return *this;
	}
	quantity& operator/=(quantity const & rhs)
	{
		value /= rhs.value;
		units /= rhs.units;
		return *this;
	}

	template<typename U>
	typename std::enable_if<std::is_convertible<U,value_type>::value, quantity&>::type operator*=(U const rhs)
	{
		value *= rhs;
		return *this;
	}
	template<typename U>
	typename std::enable_if<std::is_convertible<U,value_type>::value, quantity&>::type operator/=(U const rhs)
	{
		value /= rhs;
		return *this;
	}

	const quantity operator+(quantity const & rhs) const
	{
		return quantity(*this) += rhs;
	}
	const quantity operator-(quantity const & rhs) const
	{
		return quantity(*this) -= rhs;
	}
	const quantity operator*(quantity const & rhs) const
	{
		return quantity(*this) *= rhs;
	}
	const quantity operator/(quantity const & rhs) const
	{
		return quantity(*this) /= rhs;
	}
	template<typename U>
	typename std::enable_if<std::is_convertible<U,value_type>::value, const quantity>::type operator*(U const rhs) const
	{
		return quantity(*this) *= rhs;
	}
	template<typename U>
	typename std::enable_if<std::is_convertible<U,value_type>::value, const quantity>::type operator/(U const rhs) const
	{
		return quantity(*this) /= rhs;
	}

	bool operator==(quantity const & rhs) const
	{
		return value == rhs.value
			&& units == rhs.units;
	}
	bool operator!=(quantity const & rhs) const
	{
		return !(*this == rhs);
	}
	bool operator>(quantity const & rhs) const
	{
		if (units != rhs.units)
			throw std::runtime_error("Comparing incompatible units.");
		return value > rhs.value;
	}
	bool operator>=(quantity const & rhs) const
	{
		if (units != rhs.units)
			throw std::runtime_error("Comparing incompatible units.");
		return value >= rhs.value;
	}
	bool operator<(quantity const & rhs) const
	{
		if (units != rhs.units)
			throw std::runtime_error("Comparing incompatible units.");
		return value < rhs.value;
	}
	bool operator<=(quantity const & rhs) const
	{
		if (units != rhs.units)
			throw std::runtime_error("Comparing incompatible units.");
		return value <= rhs.value;
	}
};

template<typename T>
inline quantity<T> pow(quantity<T> q, int power)
{
	return quantity<T>(std::pow(q.value, power), pow(q.units, power));
}

template<typename T, typename U>
inline typename std::enable_if<std::is_convertible<T,U>::value, const quantity<U>>::type operator*(T const v, quantity<U> const & q)
{
	return q * v;
}
template<typename T>
inline const quantity<T> operator/(T const v, quantity<T> const & q)
{
	return quantity<T>(v / q.value, dimensions::dimensionless / q.units);
}

}} // namespace nbl::units

#endif
