#ifndef __UNIT_SYSTEM_H_
#define __UNIT_SYSTEM_H_

/*
 * A numerical choice for fundamental and derived units.
 *
 * TODO: use variable templates here, but VS does not seem to like that yet.
 */

#include "quantity.h"
#include "unit_parser.h"

namespace nbl { namespace units {

// Base units
constexpr quantity<double> dimensionless{ 1, dimensions::dimensionless };
constexpr quantity<double> eV{ 1, dimensions::energy };      // electronvolts
constexpr quantity<double> nm{ 1, dimensions::length };      // nanometers
constexpr quantity<double> ns{ 1, dimensions::time };        // nanoseconds
constexpr quantity<double> K { 1, dimensions::temperature }; // Kelvin
constexpr quantity<double> e { 1, dimensions::charge };      // electron charge


// Auxiliary units
constexpr quantity<double> s  { 1e9, dimensions::time };         // second
constexpr quantity<double> m  { 1e9, dimensions::length };       // meter
constexpr quantity<double> cm { 1e7, dimensions::length };       // centimeter
constexpr quantity<double> C  { 6.2415e18, dimensions::charge }; // Coulomb
constexpr quantity<double> g  { 6.2415e15, dimensions::mass };   // gram

inline unit_parser<double> default_unit_parser()
{
	unit_parser<double> p;

	p.add_unit("eV", units::eV);
	p.add_unit("nm", units::nm);
	p.add_unit("K",  units::K);
	p.add_unit("s",  units::s);
	p.add_unit("m",  units::m);
	p.add_unit("cm", units::cm);
	p.add_unit("g",  units::g);
	p.add_unit("C",  units::C);
	p.add_unit("radian", units::dimensionless);

	return p;
}

}} // namespace nbl::units

#endif
