#include <cstring>

namespace nbl { namespace util {

template<typename T, bool gpu_flag>
CPU table_2D<T, gpu_flag> table_2D<T, gpu_flag>::create(
	real x_min, real x_max, int width,
	real y_min, real y_max, int height,
	T* data)
{
	table_2D<T, gpu_flag> table;
	detail::table_2D_factory<T, gpu_flag>::allocate(table, width, height);
	if (data != nullptr)
		detail::table_2D_factory<T, gpu_flag>::set(table, data);

	table._x_min = x_min;
	table._x_step = (width - 1) / (x_max - x_min);
	table._y_min = y_min;
	table._y_step = (height - 1) / (y_max - y_min);

	return table;
}

template<typename T, bool gpu_flag>
template<bool source_gpu_flag>
CPU table_2D<T, gpu_flag> table_2D<T, gpu_flag>::create(
	table_2D<T, source_gpu_flag> const & source)
{
	table_2D<T, gpu_flag> table;
	detail::table_2D_factory<T, gpu_flag>::allocate(table, source._width, source._height);
	detail::table_2D_factory<T, gpu_flag>::set(table, source);
	table._x_min = source._x_min;
	table._x_step = source._x_step;
	table._y_min = source._y_min;
	table._y_step = source._y_step;
	return table;
}

template<typename T, bool gpu_flag>
CPU void table_2D<T, gpu_flag>::destroy(table_2D<T, gpu_flag> & table)
{
	detail::table_2D_factory<T, gpu_flag>::free(table);
}

template<typename T, bool gpu_flag>
template<typename callback_function>
CPU void table_2D<T, gpu_flag>::mem_scope(callback_function callback)
{
	detail::table_2D_factory<T, gpu_flag>::mem_scope(*this, callback);
}

template<typename T, bool gpu_flag>
template<bool other_gpu_flag>
CPU void table_2D<T, gpu_flag>::set(table_2D<T, other_gpu_flag> const & source)
{
	if (_width != source._width || _height != source._height)
		throw std::runtime_error("Size of target table does not match source table");
	detail::table_2D_factory<T, gpu_flag>::set(*this, source);
}

template<typename T, bool gpu_flag>
PHYSICS T & table_2D<T, gpu_flag>::operator()(int i, int j)
{
	return *(reinterpret_cast<T*>(reinterpret_cast<uint8_t*>(_data) + i * _pitch) + j);
}

template<typename T, bool gpu_flag>
PHYSICS T const & table_2D<T, gpu_flag>::operator()(int i, int j) const
{
	return *(reinterpret_cast<T*>(reinterpret_cast<uint8_t*>(_data) + i * _pitch) + j);
}

template<typename T, bool gpu_flag>
PHYSICS T table_2D<T, gpu_flag>::get(real x, real y) const
{
	const real x_index = (x - _x_min) * _x_step;
	const real y_index = (y - _y_min) * _y_step;

	const int low_x = clampi(floorr(x_index), 0, _width - 2);
	const int low_y = clampi(floorr(y_index), 0, _height - 2);

	const T v00 = (*this)(low_x, low_y);
	const T v10 = (*this)(low_x + 1, low_y);
	const T v01 = (*this)(low_x, low_y + 1);
	const T v11 = (*this)(low_x + 1, low_y + 1);

	const real frac_x = x_index - low_x;
	const real frac_y = y_index - low_y;

	return (1 - frac_x)*(1 - frac_y)*v00
		+ frac_x*(1 - frac_y)*v10
		+ (1 - frac_x)*frac_y*v01
		+ frac_x*frac_y*v11;
}

template<typename T, bool gpu_flag>
PHYSICS T table_2D<T, gpu_flag>::get_rounddown(real x, real y) const
{
	const real x_index = (x - _x_min) * _x_step;
	const real y_index = (y - _y_min) * _y_step;
	const int low_x = clampi(floorr(x_index), 0, _width - 1);
	const int low_y = clampi(floorr(y_index), 0, _height - 1);
	return (*this)(low_x, low_y);
}

template<typename T, bool gpu_flag>
PHYSICS T table_2D<T, gpu_flag>::get_nearest(real x, real y) const
{
	const real x_index = (x - _x_min) * _x_step;
	const real y_index = (y - _y_min) * _y_step;
	const int near_x = clampi(rintr(x_index), 0, _width - 1);
	const int near_y = clampi(rintr(y_index), 0, _height - 1);
	return (*this)(near_x, near_y);
}

template<typename T, bool gpu_flag>
PHYSICS int table_2D<T, gpu_flag>::width() const
{
	return _width;
}

template<typename T, bool gpu_flag>
PHYSICS int table_2D<T, gpu_flag>::height() const
{
	return _height;
}

template<typename T, bool gpu_flag>
CPU void table_2D<T, gpu_flag>::set_scale(
      real x_min, real x_max,
      real y_min, real y_max)
{
	_x_min = x_min;
	_x_step = (_width - 1) / (x_max - x_min);
	_y_min = y_min;
	_y_step = (_height - 1) / (y_max - y_min);
}

template<typename T, bool gpu_flag>
template<bool gpu2>
CPU typename std::enable_if<!gpu2, T*>::type table_2D<T, gpu_flag>::data()
{
	return _data;
}

namespace detail
{
	template<typename T>
	struct table_2D_factory<T, false>
	{
		inline static CPU void allocate(table_2D<T, false> & table, int width, int height)
		{
			table._width = width;
			table._height = height;
			table._data = new T[width * height];
			table._pitch = height * sizeof(T);
		}

		inline static CPU void set(table_2D<T, false> & table, T* data)
		{
			memcpy(table._data, data,
				table._width * table._height * sizeof(T));
		}
		inline static CPU void set(table_2D<T, false> & target, table_2D<T, false> const & source)
		{
			memcpy(target._data, source._data,
				target._width * target._height * sizeof(T));
		}
#if CUDA_COMPILER_AVAILABLE
		inline static CPU void set(table_2D<T, false> & target, table_2D<T, true> const & source)
		{
			cudaMemcpy2D(target._data, target._pitch,
				source._data, source._pitch,
				target._height*sizeof(T), target._width,
				cudaMemcpyDeviceToHost);
		}
#endif // CUDA_COMPILER_AVAILABLE

		template<typename callback_function>
		inline static CPU void mem_scope(table_2D<T, false> & table, callback_function callback)
		{
			// Make indirect array
			T** host_pp = new T*[table._width];
			for (int x = 0; x < table._width; ++x)
				host_pp[x] = reinterpret_cast<T*>(reinterpret_cast<uint8_t*>(table._data) + x * table._pitch);

			callback(host_pp);

			delete[] host_pp;
		}

		inline static CPU void free(table_2D<T, false> & table)
		{
			delete[] table._data;
			table._data = nullptr;
			table._width = 0;
			table._height = 0;
			table._pitch = 0;
		}
	};

#if CUDA_COMPILER_AVAILABLE
	template<typename T>
	struct table_2D_factory<T, true>
	{
		inline static CPU void allocate(table_2D<T, true> & table, int width, int height)
		{
			table._width = width;
			table._height = height;
			cuda::cuda_new_2D(&table._data, &table._pitch, width, height);
		}

		inline static CPU void set(table_2D<T, true> & table, T* data)
		{
			cudaMemcpy2D(table._data, table._pitch,
				data, table._height*sizeof(T),
				table._height*sizeof(T), table._width,
				cudaMemcpyHostToDevice);
		}
		inline static CPU void set(table_2D<T, true> & target, table_2D<T, false> const & source)
		{
			cudaMemcpy2D(target._data, target._pitch,
				source._data, source._pitch,
				target._height*sizeof(T), target._width,
				cudaMemcpyHostToDevice);
		}
		inline static CPU void set(table_2D<T, true> & target, table_2D<T, true> const & source)
		{
			// We only support CPU-GPU clones for now.
			// To go the other way, the user expects us to do a peer-to-peer copy,
			// which we can't do because we don't want to rely on unified addressing
			// and this class does not keep track of the GPU its data resides on.
			static_assert(true, "GPU-GPU clones not supported yet.");
		}

		template<typename callback_function>
		inline static CPU void mem_scope(table_2D<T, true> & table, callback_function callback)
		{
			cuda::cuda_mem_scope_2D<T>(table._data, table._pitch, table._width, table._height, callback);
		}

		inline static CPU void free(table_2D<T, true> & table)
		{
			cudaFree(table._data);
			table._data = nullptr;
			table._width = 0;
			table._height = 0;
			table._pitch = 0;
		}
	};
#endif // CUDA_COMPILER_AVAILABLE
} // namespace detail

}} // namespace nbl::util
