#ifndef __DIMENSION_H_
#define __DIMENSION_H_

/*
 * We have chosen some fundamental unit dimensions, from which others are derived
 * as a combination of powers of fundamental units. The dimension struct keeps track
 * of these powers, only integer powers for now.
 * This struct keeps track of dimensions, such as "energy" and "length".
 * What consitites "one energy" (e.g. electronvolt, joule) is defined elsewhere.
 *
 * The bottom of this file defines some standard dimensions to be used elsewhere.
 */

namespace nbl { namespace units {

struct dimension
{
	int energy;
	int length;
	int time;
	int temperature;
	int charge;

	dimension& operator*=(dimension const & rhs)
	{
		energy += rhs.energy;
		length += rhs.length;
		time += rhs.time;
		temperature += rhs.temperature;
		charge += rhs.charge;
		return *this;
	}
	dimension& operator/=(dimension const & rhs)
	{
		energy -= rhs.energy;
		length -= rhs.length;
		time -= rhs.time;
		temperature -= rhs.temperature;
		charge -= rhs.charge;
		return *this;
	}

	const dimension operator*(dimension const & rhs) const
	{
		return dimension(*this) *= rhs;
	}
	const dimension operator/(dimension const & rhs) const
	{
		return dimension(*this) /= rhs;
	}

	bool operator==(dimension const & rhs) const
	{
		return energy == rhs.energy
			&& length == rhs.length
			&& time == rhs.time
			&& temperature == rhs.temperature
			&& charge == rhs.charge;
	}
	bool operator!=(dimension const & rhs) const
	{
		return !(*this == rhs);
	}
};

inline dimension pow(dimension dim, int power)
{
	return dimension
	{
		dim.energy * power,
		dim.length * power,
		dim.time * power,
		dim.temperature * power,
		dim.charge * power
	};
}

namespace dimensions
{
	// Fundamental dimensions
	constexpr dimension dimensionless{ 0, 0, 0, 0, 0 };
	constexpr dimension energy       { 1, 0, 0, 0, 0 };
	constexpr dimension length       { 0, 1, 0, 0, 0 };
	constexpr dimension time         { 0, 0, 1, 0, 0 };
	constexpr dimension temperature  { 0, 0, 0, 1, 0 };
	constexpr dimension charge       { 0, 0, 0, 0, 1 };

	// Derived dimensions
	constexpr dimension mass         { 1, -2, 2, 0, 0 };
	constexpr dimension acceleration { 0, 1, -2, 0, 0 };
	constexpr dimension force        { 1, -1, 0, 0, 0 };
	constexpr dimension area         { 0, 2, 0, 0, 0 };
	constexpr dimension volume       { 0, 3, 0, 0, 0 };
	constexpr dimension n_density    { 0, -3, 0, 0, 0 }; // Number density
	constexpr dimension m_density    { 1, -5, 2, 0, 0 }; // Mass density
}

}} // namespace nbl::units

#endif
