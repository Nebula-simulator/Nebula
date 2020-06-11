#ifndef __TIME_LOG_H_
#define __TIME_LOG_H_

#include <chrono>
#include <iomanip>

/**
 * \brief Very simple class to print some time logging data after a simulation
 *        finishes.
 */

struct time_log
{
public:
	/// Start timer
	void start()
	{
		last_start = std::chrono::steady_clock::now();
	}

	/**
	 * \brief Stop timer
	 *
	 * The time elapsed since the last call to start() is added to the list.
	 *
	 * \param name Name to use for printing
	 */
	void stop(std::string const & name)
	{
		auto end = std::chrono::steady_clock::now();
		data.push_back({name, end - last_start});
	}

	/**
	 * \brief Print a text summary
	 *
	 * \param out Stream to print the summary to.
	 */
	void print(std::ostream& out)
	{
		// Measure the length of the first column
		const size_t len_col1 = std::max_element(data.begin(), data.end(),
			[](data_t const & d1, data_t const & d2) -> bool
			{ return d1.first.size() < d2.first.size(); })->first.size();

		// Print the data
		out << "Timing (seconds):\n";
		for (auto&& d : data)
		{
			const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(d.second);
			out << std::setw(len_col1) << d.first << ": "
				<< std::fixed << std::setprecision(3) << (ms.count()/1000.)
				<< '\n';
		}
	}

private:
	using data_t = std::pair<std::string, std::chrono::steady_clock::duration>;

	std::chrono::steady_clock::time_point last_start;
	std::vector<data_t> data;
};

#endif // __TIME_LOG_H_
