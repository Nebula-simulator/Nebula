#ifndef __UNIT_PARSER_H_
#define __UNIT_PARSER_H_

/*
 * Class to parse number - unit strings.
 */

#include <string>
#include <map>
#include <cctype>

namespace nbl { namespace units {

template<typename T>
class unit_parser
{
public:
	using value_type = quantity<T>;

	value_type parse_unit(std::string const & text) const
	{
		return parse_unit(text.cbegin(), text.cend());
	}

	value_type parse_value_unit(std::string const & text) const
	{
		// Parse number, simply by throwing the whole string into std::stod.
		std::size_t next_idx;
		const double value = std::stod(text, &next_idx);

		// The part that std::stod could not parse, is the unit.
		return value * parse_unit(text.cbegin() + next_idx, text.cend());
	}

	void add_unit(std::string const & unit_name, value_type unit_value)
	{
		unit_map.insert({ unit_name, unit_value });
	}

private:
	std::map<std::string, value_type> unit_map;

	/*
	 * This function parses a unit in the form of "g/cm^3", possibly with a leading * or /.
	 */
	template<typename iterator_type>
	value_type parse_unit(iterator_type current_position, iterator_type end) const
	{
		std::string buffer;
		value_type final_unit{ 1, dimensions::dimensionless };

		while (current_position != end)
		{
			// Get operator (* or /) -- overall multiplier for the power.
			int power = (*current_position == '/' ? -1 : 1);
			if (is_new_unit(*current_position))
				++current_position;


			// Get unit string (e.g. 'eV', 'nm', ...)
			buffer.clear();
			for (; current_position != end; ++current_position)
			{
				if (is_power(*current_position) || is_new_unit(*current_position))
					break;
				buffer.push_back(*current_position);
			}
			auto this_unit_iterator = unit_map.find(buffer);
			if (this_unit_iterator == unit_map.end())
				throw std::runtime_error("Unknown unit '" + buffer + "'.");
			const auto this_unit = this_unit_iterator->second;


			// Get trailing power, if any
			if (current_position != end && is_power(*current_position))
			{
				buffer.clear();
				for (++current_position; current_position != end; ++current_position)
				{
					if (is_new_unit(*current_position))
						break;
					buffer.push_back(*current_position);
				}
				power *= std::stoi(buffer);
			}


			final_unit *= pow(this_unit, power);
			// Note that, at this point, is_new_unit(*current_position) == true or current_position == end.
		}
		return final_unit;
	}

	static bool is_new_unit(char c)
	{
		return c == '*' || c == '/';
	}
	static bool is_power(char c)
	{
		return c == '^';
	}
};

}} // namespace nbl::units

#endif
