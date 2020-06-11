#ifndef __CLI_PARAMS_H_
#define __CLI_PARAMS_H_

/**
 * \brief Very simple command-line parameter parser.
 *
 * Types of arguments supported are:
 *  * `--key=value`
 *  * `--key value`
 *  * Anything that doesn't match this pattern are "positional arguments".
 */

#include <string>
#include <vector>
#include <map>
#include <sstream>

class cli_params
{
public:
	/**
	 * \brief Constructor.
	 *
	 * \param param_string String to be printed for the documentation listing.
	 */
	cli_params(std::string const & param_string = "[options]")
		: param_string(param_string)
	{}

	/**
	 * \brief Add key-value option pair.
	 */
	template<typename T>
	void add_option(
		std::string const & key,
		std::string const & description,
		T const & default_value)
	{
		std::stringstream def;
		def << default_value;
		options[key] = { description, def.str() };
	}

	/**
	 * \brief Parse command-line options.
	 *
	 * Parameters are the same as what you get in the C++ `main()` function.
	 *
	 * This function parses the arguments and stores them in private variables.
	 *
	 * \param argc Number of command-line arguments.
	 * \param argv The command-line arguments themselves.
	 */
	void parse(int argc, char**argv)
	{
		program_name = argv[0];
		for (int i = 1; i < argc; ++i)
		{
			std::string item = argv[i];
			if (item[0] == '-' && item[1] == '-')
			{
				std::string key;
				std::string value;

				auto split = item.find("=", 2);
				if (split == std::string::npos)
				{
					key = item.substr(2, item.size()-2);
					value = (++i < argc ? argv[i] : "");
					if (value == "=")
						value = (++i < argc ? argv[i] : "");
				}
				else
				{
					key = item.substr(2, split-2);
					value = item.substr(split+1);
				}

				auto map_it = options.find(key);
				if (map_it == options.end())
					throw std::runtime_error("Unexpected key '" + key + "'");
				map_it->second.value = value;
			}
			else
			{
				positional.push_back(item);
			}
		}
	}

	/**
	 * \brief Print documentation string
	 */
	void print_usage(std::ostream & out)
	{
		out << "Usage: " << program_name << " " << param_string << "\n"
			<< "Options:\n";
		for (auto const & opt : options)
		{
			out << "\t--" << opt.first << ": " << opt.second.description
				<< " [" << opt.second.value << "]\n";
		}
	}

	/**
	 * \brief Get the program name, i.e. argv[0].
	 */
	std::string const & get_program_name() const
	{
		return program_name;
	}

	/**
	 * \brief Get a flag.
	 *
	 * \param key Flag name
	 */
	template<typename T>
	T get_flag(std::string const & key) const
	{
		return options.at(key).get_as<T>();
	}

	/**
	 * \brief Get the positional arguments (in the correct order).
	 *
	 * The "positional" arguments are arguments that don't follow the --key=value format.
	 */
	std::vector<std::string> const & get_positional() const
	{
		return positional;
	}

private:
	struct option
	{
		std::string description;
		std::string value;

		template<typename T>
		T get_as() const
		{
			T ret;
			std::stringstream(value) >> ret;
			return ret;
		}
	};

	std::string program_name;
	std::string param_string;
	std::map<std::string, option> options;
	std::vector<std::string> positional;
};

#endif // __CLI_PARAMS_H_
