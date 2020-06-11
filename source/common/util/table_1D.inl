#include <cstring>

namespace nbl { namespace util {

template<typename T, bool gpu_flag>
CPU table_1D<T, gpu_flag> table_1D<T, gpu_flag>::create(
	real x_min, real x_max, int n,
	T* data)
{
	table_1D<T, gpu_flag> table;
	detail::table_1D_factory<T, gpu_flag>::allocate(table, n);
	if (data != nullptr)
		detail::table_1D_factory<T, gpu_flag>::set(table, data);
	table._x_min = x_min;
	table._x_step = (n - 1) / (x_max - x_min);
	return table;
}

template<typename T, bool gpu_flag>
template<bool source_gpu_flag>
CPU table_1D<T, gpu_flag> table_1D<T, gpu_flag>::create(
	table_1D<T, source_gpu_flag> const & source)
{
	table_1D<T, gpu_flag> table;
	detail::table_1D_factory<T, gpu_flag>::allocate(table, source._n);
	detail::table_1D_factory<T, gpu_flag>::set(table, source);
	table._x_min = source._x_min;
	table._x_step = source._x_step;
	return table;
}

template<typename T, bool gpu_flag>
CPU void table_1D<T, gpu_flag>::destroy(table_1D<T, gpu_flag> & table)
{
	detail::table_1D_factory<T, gpu_flag>::free(table);
}

template<typename T, bool gpu_flag>
template<typename callback_function>
CPU void table_1D<T, gpu_flag>::mem_scope(callback_function callback)
{
	detail::table_1D_factory<T, gpu_flag>::mem_scope(*this, callback);
}

template<typename T, bool gpu_flag>
template<bool other_gpu_flag>
CPU void table_1D<T, gpu_flag>::set(table_1D<T, other_gpu_flag> const & source)
{
	if (_n != source._n)
		throw std::runtime_error("Size of target table does not match source table");
	detail::table_1D_factory<T, gpu_flag>::set(*this, source);
}

template<typename T, bool gpu_flag>
PHYSICS T & table_1D<T, gpu_flag>::operator()(int i)
{
	return _data[i];
}
template<typename T, bool gpu_flag>
PHYSICS T const & table_1D<T, gpu_flag>::operator()(int i) const
{
	return _data[i];
}

template<typename T, bool gpu_flag>
PHYSICS T table_1D<T, gpu_flag>::get(real x) const
{
	const real x_index = (x - _x_min) * _x_step;

	const int low_index = clampi(floorr(x_index), 0, _n - 2);
	const T low_value = _data[low_index];
	const T high_value = _data[low_index + 1];

	const real frac_index = x_index - low_index;
	return (1 - frac_index)*low_value + frac_index*high_value;

/*
Original implementation does not handle infinities correctly.
	return low_value + (x_index - low_index) * (high_value - low_value);
*/
}

template<typename T, bool gpu_flag>
PHYSICS int table_1D<T, gpu_flag>::width() const
{
	return _n;
}

template<typename T, bool gpu_flag>
CPU void table_1D<T, gpu_flag>::set_scale(real x_min, real x_max)
{
	_x_min = x_min;
	_x_step = (_n - 1) / (x_max - x_min);
}

template<typename T, bool gpu_flag>
PHYSICS real table_1D<T, gpu_flag>::get_scalemin() const
{
	return _x_min;
}

template<typename T, bool gpu_flag>
PHYSICS real table_1D<T, gpu_flag>::get_scalemax() const
{
	return _x_min + (_n - 1) / _x_step;
}

template<typename T, bool gpu_flag>
template<bool gpu2>
CPU typename std::enable_if<!gpu2, T*>::type table_1D<T, gpu_flag>::data()
{
	return _data;
}

namespace detail
{
	template<typename T>
	struct table_1D_factory<T, false>
	{
		inline static CPU void allocate(table_1D<T, false> & table, int n)
		{
			table._n = n;
			table._data = new T[n];
		}

		inline static CPU void set(table_1D<T, false> & table, T* data)
		{
			memcpy(table._data, data, table._n * sizeof(T));
		}
		inline static CPU void set(table_1D<T, false> & target, table_1D<T, false> const & source)
		{
			memcpy(target._data, source._data, target._n*sizeof(T));
		}
#if CUDA_COMPILER_AVAILABLE
		inline static CPU void set(table_1D<T, false> & target, table_1D<T, true> const & source)
		{
			cudaMemcpy(target._data, source._data, target._n*sizeof(T), cudaMemcpyDeviceToHost);
		}
#endif // CUDA_COMPILER_AVAILABLE

		template<typename callback_function>
		inline static CPU void mem_scope(table_1D<T, false> & table, callback_function callback)
		{
			callback(table._data);
		}

		inline static CPU void free(table_1D<T, false> & table)
		{
			delete[] table._data;
			table._data = nullptr;
			table._n = 0;
		}
	};

#if CUDA_COMPILER_AVAILABLE
	template<typename T>
	struct table_1D_factory<T, true>
	{
		inline static CPU void allocate(table_1D<T, true> & table, int n)
		{
			table._n = n;
			cuda::cuda_new<T>(&table._data, n);
		}

		inline static CPU void set(table_1D<T, true> & table, T* data)
		{
			cudaMemcpy(table._data, data, table._n*sizeof(T), cudaMemcpyHostToDevice);
		}
		inline static CPU void set(table_1D<T, true> & target, table_1D<T, false> const & source)
		{
			cudaMemcpy(target._data, source._data, target._n*sizeof(T), cudaMemcpyHostToDevice);
		}
		inline static CPU void set(table_1D<T, true> & target, table_1D<T, true> const & source)
		{
			// We only support CPU-GPU clones for now.
			// To go the other way, the user expects us to do a peer-to-peer copy,
			// which we can't do because we don't want to rely on unified addressing
			// and this class does not keep track of the GPU its data resides on.
			static_assert(true, "GPU-GPU clones not supported yet.");
		}

		template<typename callback_function>
		inline static CPU void mem_scope(table_1D<T, true> & table, callback_function callback)
		{
			cuda::cuda_mem_scope<T>(table._data, table._n, callback);
		}

		inline static CPU void free(table_1D<T, true> & table)
		{
			cudaFree(table._data);
			table._data = nullptr;
			table._n = 0;
		}
	};
#endif // CUDA_COMPILER_AVAILABLE
} // namespace detail

}} // namespace nbl::util
