#ifndef __TABLE_3D_H_
#define __TABLE_3D_H_

namespace nbl { namespace util {

namespace detail
{
	/**
	 * \brief Responsible for managing the memory on the CPU or GPU device.
	 *
	 * Reason for putting this in a separate struct is that we can't
	 * partially specialize member functions of the table_3D class.
	 */
	template<typename T, bool gpu_flag>
	struct table_3D_factory;
}

/**
 * \brief 3D table to be used in simulation code.
 *
 * It doesn't just store data, but also information about the "axes".
 *
 * The table is allocated by the static {@link create} function and deallocated
 * by the {@link destroy} function. It is up to the user to call these functions
 * correctly! The reason is that these tables are expected to be passed to the
 * GPU as parameters, so we don't want to automatically call constructors and
 * destructors.
 *
 * \tparam T        Datatype stored (usually ::real)
 * \tparam gpu_flag Set to true if data should be stored on the GPU, false for CPU.
 */
template<typename T, bool gpu_flag>
class table_3D
{
public:
	/**
	 * \brief Create a table, allocating memory for the data.
	 *
	 * Throws `std::bad_alloc` exception if memory allocation fails.
	 *
	 * \param x_min  Lower value of the physical range along the first dimension
	 * \param x_max  Upper value of the physical range along the first dimension
	 * \param width  Number of data points in x direction
	 * \param y_min  Lower value of the physical range along the second dimension
	 * \param y_max  Upper value of the physical range along the second dimension
	 * \param height Number of data points in y direction
	 * \param z_min  Lower value of the physical range along the third dimension
	 * \param z_max  Upper value of the physical range along the third dimension
	 * \param depth  Number of data points in z direction
	 * \param data   Optionally: provide data. Order is such that
	 *               `data[x,y,z] = data + z*width*height + y*width + x`.
	 */
	static CPU table_3D<T, gpu_flag> create(
		real x_min, real x_max, int width,
		real y_min, real y_max, int height,
		real z_min, real z_max, int depth,
		T* data = nullptr);

	/**
	 * \brief Create a table, as a deep copy of another one.
	 */
	template<bool source_gpu_flag>
	static CPU table_3D<T, gpu_flag> create(table_3D<T, source_gpu_flag> const & source);

	/**
	 * \brief Destroy a table, deallocating the data.
	 */
	static CPU void destroy(table_3D<T, gpu_flag> & arr);

	/**
	 * \brief Set or get data using a callback function.
	 *
	 * Similar to ::nbl::cuda::cuda_mem_scope_3D. If applicable, it copies the
	 * data from the GPU to the CPU. It calls the \p callback function with a
	 * CPU-accessible pointer to this data. Finally, if applicable, it copies
	 * the data back to the GPU.
	 *
	 * The indexing convention for the callback function is `array[z][y][x]`.
	 *
	 * \param callback Callback function. Should have signature `void(T***)`.
	 */
	template<typename callback_function>
	CPU void mem_scope(callback_function callback);

	/**
	 * \brief Copy the data from another table to this one.
	 */
	template<bool other_gpu_flag>
	CPU void set(table_3D<T, other_gpu_flag> const & source);

	/**
	 * \brief Direct read-write access to data. No bounds checking is done.
	 *
	 * \param i Index to be accessed in the first dimension.
	 * \param j Index to be accessed in the second dimension.
	 * \param k Index to be accessed in the third dimension.
	 */
	inline PHYSICS T & operator()(int i, int j, int k);
	/**
	 * \brief Direct read-only access to data. No bounds checking is done.
	 *
	 * \param i Index to be accessed in the first dimension.
	 * \param j Index to be accessed in the second dimension.
	 * \param k Index to be accessed in the third dimension.
	 */
	inline PHYSICS T const & operator()(int i, int j, int k) const;

	/**
	 * \brief Get value at some (x,y,z) coordinate, with linear interpolation
	 * and extrapolation.
	 *
	 * If \p x is between `x_min` and `x_max` provided to the constructor, this
	 * function performs linear interpolation. If it is outside that range, this
	 * function performs linear extrapolation. Same for \p y and \p z.
	 *
	 * \param x Position along the first dimension to get the data for.
	 * \param y Position along the second dimension to get the data for.
	 * \param z Position along the third dimension to get the data for.
	 */
	inline PHYSICS T get(real x, real y, real z) const;
	/**
	 * \brief Get value at some (x,y,z) coordinate, without interpolation.
	 *
	 * Converts position to indices, rounds them down, and returns the
	 * non-interpolated value. If \p x, \p y or \p z are below the lowest value,
	 * rounds up to index 0.
	 *
	 * \param x Position along the first dimension to get the data for.
	 * \param y Position along the second dimension to get the data for.
	 * \param z Position along the third dimension to get the data for.
	 */
	inline PHYSICS T get_rounddown(real x, real y, real z) const;
	/**
	 * \brief Get value at some (x,y,z) coordinate, without interpolation.
	 *
	 * Converts position to indices, rounds them to nearest, and returns the
	 * non-interpolated value.
	 *
	 * \param x Position along the first dimension to get the data for.
	 * \param y Position along the second dimension to get the data for.
	 * \param z Position along the third dimension to get the data for.
	 */
	inline PHYSICS T get_nearest(real x, real y, real z) const;

	/**
	 * \brief Get the width of this table.
	 */
	inline PHYSICS int width() const;

	/**
	 * \brief Get the height of this table.
	 */
	inline PHYSICS int height() const;

	/**
	 * \brief Get the depth of this table.
	 */
	inline PHYSICS int depth() const;

	/**
	 * \brief Set the scale.
	 *
	 * The data itself is preserved.
	 */
	inline CPU void set_scale(
		   real x_min, real x_max,
		   real y_min, real y_max,
		   real z_min, real z_max);

	/**
	 * \brief Get write access to the underlying data array.
	 *
	 * This is only possible if \p gpu_flag is false.
	 * The underlying data array is stored in C order, that is, the fastest
	 * changing index is the last one (running from 0 to depth).
	 */
	template<bool gpu2=gpu_flag>
	inline CPU typename std::enable_if<!gpu2, T*>::type data();

private:
	T* _data = nullptr;
	size_t _pitch = 0;

	int _width = 0;
	real _x_min = 0;
	real _x_step = 0;

	int _height = 0;
	real _y_min = 0;
	real _y_step = 0;

	int _depth = 0;
	real _z_min = 0;
	real _z_step = 0;

	friend struct detail::table_3D_factory<T, gpu_flag>;
	friend struct detail::table_3D_factory<T, !gpu_flag>;
	friend class table_3D<T, !gpu_flag>;
};

}} // namespace nbl::util

#include "table_3D.inl"

#endif // __TABLE_3D_H_
