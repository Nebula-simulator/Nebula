#ifndef __SCATTER_LIST_H_
#define __SCATTER_LIST_H_


#include "../common/tuple.h"
#include "events.h"
#include "particle.h"

// TODO namespace

/**
 * \brief Represents a list of scattering processes that can take place in a
 * material.
 *
 * \tparam scatter_types... Types of scattering events.
 *
 * Each \p scatter_type is a class that must have public member functions with
 * the following signatures:
 *   * `PHYSICS real sample_path(particle const &, random_generator&) const;`
 *
 *     Returns a free path length until the next scattering event for the given particle.
 *   * `PHYSICS void execute_event(particle_manager, particle_index, random_generator&) const;`
 *
 *     Performs the event, if chosen.
 */
template<typename... scatter_types>
class scatter_list : public nbl::tuple::tuple<scatter_types...>
{
public:
	/**
	 * \brief Default constructor
	 */
	scatter_list() = default;

	/**
	 * \brief Constructor
	 */
	CPU scatter_list(scatter_types... sc)
		: nbl::tuple::tuple<scatter_types...>(sc...)
	{}

	/**
	 * \brief Get the scattering type at the given index
	 */
	template <size_t I>
	using type_at_index = typename nbl::type_at_index<I, scatter_types...>;

	/**
	 * \brief Return the number of scattering event types
	 */
	static constexpr size_t size() { return sizeof...(scatter_types); }

	/**
	 * \brief Print information summary
	 */
	static void print_info(std::ostream& stream)
	{
		(void)std::initializer_list<int>{(scatter_types::print_info(stream), 0)... };
	}

	/**
	 * \brief Get maximal energy, in eV
	 */
	real get_max_energy() const
	{
		energy_visitor visitor;
		nbl::tuple::for_each(*this, visitor);
		return visitor.max_energy;
	}

/*
	// Generic lambdas... would be nice, right?
	// But we chose not to rely on C++14 yet...

	template<typename rng_t>
	inline PHYSICS scatter_event sample_path(particle const & this_particle, rng_t & rng) const
	{

		scatter_event evt{ 0, std::numeric_limits<real>::infinity() };
		nbl::tuple::for_each(*this, [&](auto scatter, size_t i)
			{
				auto path = scatter.sample_path(this_particle, rng);
				if (path < evt.distance)
					evt = { uint8_t(i+1), path };
			});
		return evt;
	}

	template<typename particle_manager, typename rng_t>
	inline PHYSICS void execute(uint8_t i,
		particle_manager& particle_mgr,
		typename particle_manager::particle_index_t particle_idx,
		rng_t& rng) const
	{
		nbl::tuple::visit_at(*this, i-1,
			[&](auto scatter)
			{ scatter.execute(particle_mgr, particle_idx, rng); });
	}
*/

	// C++11 workarounds
	/**
	 * \brief Returns a ::scatter_event describing which event takes place
	 *        first, and at what distance
	 */
	template<typename rng_t>
	inline PHYSICS scatter_event sample_path(particle const & this_particle, rng_t & rng) const
	{
		sample_visitor<rng_t> visitor{ this_particle, rng, { 0, std::numeric_limits<real>::infinity() }};
		nbl::tuple::for_each(*this, visitor);
		return visitor.evt;
	}

	/**
	 * \brief Executes the desired event.
	 *
	 * \param i            Index of the scattering event type (typically, the
	 *                     index returned by #sample_path).
	 * \param particle_mgr Particle manager used by the simulation.
	 * \param particle_idx Index of the particle for which the event is to be
	 *                     executed.
	 * \param rng          A random number generator.
	 */
	template<typename particle_manager, typename rng_t>
	inline PHYSICS void execute(uint8_t i,
		particle_manager& particle_mgr,
		typename particle_manager::particle_index_t particle_idx,
		rng_t& rng) const
	{
		execute_visitor<particle_manager, rng_t> visitor{ particle_mgr, particle_idx, rng };
		nbl::tuple::visit_at(*this, i-1, visitor);
	}

private:
	struct energy_visitor
	{
		real max_energy = std::numeric_limits<real>::infinity();

		template<typename scatter_type>
		inline PHYSICS void operator()(scatter_type& scatter, size_t)
		{
			const real e = scatter.get_max_energy();
			if (e < max_energy)
				max_energy = e;
		}
	};

	template<typename rng_t>
	struct sample_visitor
	{
		particle const & this_particle;
		rng_t & rng;
		scatter_event evt;

		template<typename scatter_type>
		inline PHYSICS void operator()(scatter_type& scatter, size_t i)
		{
			auto path = scatter.sample_path(this_particle, rng);
			if (path < evt.distance)
				evt = { uint8_t(i+1), path };
		}
	};
	template<typename particle_manager, typename rng_t>
	struct execute_visitor
	{
		particle_manager& particle_mgr;
		typename particle_manager::particle_index_t& particle_idx;
		rng_t& rng;

		template<typename scatter_type>
		inline PHYSICS void operator()(scatter_type& scatter)
		{
			scatter.execute(particle_mgr, particle_idx, rng);
		}
	};
};

#endif // __SCATTER_LIST_H_
