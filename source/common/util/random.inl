#include <random>

namespace nbl { namespace util {

namespace detail
{
	template<>
	class random_state<false>
	{
	public:
		using seed_type = typename std::mt19937::result_type;
		static constexpr seed_type default_seed = std::mt19937::default_seed;

		CPU random_state(seed_type seed = default_seed)
			: _generator(seed)
		{}

		PHYSICS real unit()
		{
#if CUDA_COMPILING
			// TODO: proper error message.
			return 0;
#else
			return std::generate_canonical<real, std::numeric_limits<real>::digits>(_generator);
#endif
		}

	private:
		std::mt19937 _generator;
	};

#if CUDA_COMPILER_AVAILABLE
	#include <curand_kernel.h>

	template<>
	class random_state<true>
	{
	public:
		using seed_type = unsigned long long;
		static constexpr seed_type default_seed = 0;

		__device__ random_state(seed_type seed, seed_type sequence)
		{
			curand_init(seed, sequence, 0, &_rand_state);
		}
		PHYSICS real unit()
		{
#if CUDA_COMPILING
#if USE_DOUBLE
			return curand_uniform_double(&_rand_state);
#else // USE_DOUBLE
			return curand_uniform(&_rand_state);
#endif // USE_DOUBLE
#else // CUDA_COMPILING
			// TODO: proper error message.
			return 0;
#endif // CUDA_COMPILING
		}

	private:
		curandState _rand_state;
	};

#endif // CUDA_COMPILER_AVAILABLE
} // namespace detail

}} // namespace nbl::util
