#ifndef __CUDA_MEM_SCOPE_H_
#define __CUDA_MEM_SCOPE_H_

/**
 * \file common/cuda/cuda_mem_scope.h
 * \brief Utility functions for accessing GPU memory from CPU code
 */

namespace nbl { namespace cuda {

/**
 * \brief Utility function to read/write data in GPU memory from CPU code.
 *
 * Actions: copy GPU memory to host, call callback function to set new data,
 * copy to device. The GPU memory is an array of \p N elements of type \p T.
 *
 * \param dev_p    Pointer to GPU memory
 * \param N        Number of elements in GPU memory
 * \param callback Callback function to set the data. Should have signature `void(T*)`.
 */
template<typename T, typename callback_func>
CPU void cuda_mem_scope(T* dev_p, size_t N, callback_func callback)
{
	if (dev_p == nullptr || N == 0)
		return;

	// Copy to host.
	T* host_p;
	cudaMallocHost(&host_p, N * sizeof(T));
	cudaMemcpy(host_p, dev_p, N * sizeof(T), cudaMemcpyDeviceToHost);

	// Callback
	callback(host_p);

	// Copy back to device
	cudaMemcpy(dev_p, host_p, N * sizeof(T), cudaMemcpyHostToDevice);
	cudaFreeHost(host_p);
}

/**
 * \brief Read/write 2D array in GPU memory.
 *
 * Similar to ::cuda_mem_scope, only with 2D access.
 * The callback function should have signature `void(T** arr)`; the indexing
 * convention is `arr[x][y]`.
 *
 * \param dev_p    Pointer to GPU memory
 * \param pitch    Pitch, in bytes
 * \param width    Width (x size) of the 2D array
 * \param height   Height (y size) of the 2D array
 * \param callback The callback function. Signature should be `void(T** arr)`,
 *                 the indexing convention is `arr[x][y]`.
 */
template<typename T, typename callback_func>
CPU void cuda_mem_scope_2D(T* dev_p, size_t pitch, size_t width, size_t height, callback_func callback)
{
	if (dev_p == nullptr || pitch == 0 || width == 0 || height == 0)
		return;

	// Copy to host
	T* host_p;
	cudaMallocHost(&host_p, pitch*width);
	cudaMemcpy(host_p, dev_p, pitch*width, cudaMemcpyDeviceToHost);

	// Make indirect array
	T** host_pp = new T*[width];
	for (size_t x = 0; x < width; ++x)
		host_pp[x] = reinterpret_cast<T*>(reinterpret_cast<uint8_t*>(host_p) + x*pitch);

	// Callback
	callback(host_pp);

	// Copy back
	cudaMemcpy(dev_p, host_p, pitch*width, cudaMemcpyHostToDevice);

	delete[] host_pp;
	cudaFreeHost(host_p);
}

/**
 * \brief Read/write 3D array in GPU memory.
 *
 * Similar to ::cuda_mem_scope, only with 3D access.
 * The callback function should have signature `void(T*** arr)`; the indexing
 * convention is `arr[x][y][z]`.
 *
 * \param dev_p    Pointer to GPU memory
 * \param pitch    Pitch, in bytes
 * \param width    Width (x size) of the 3D array
 * \param height   Height (y size) of the 3D array
 * \param depth    Depth (z size) of the 3D array
 * \param callback The callback function. Signature should be `void(T*** arr)`,
 *                 the indexing convention is `arr[x][y][z]`.
 */
template<typename T, typename callback_func>
CPU void cuda_mem_scope_3D(T* dev_p, size_t pitch, size_t width, size_t height, size_t depth, callback_func callback)
{
	if (dev_p == nullptr || pitch == 0 || width == 0 || height == 0 || depth == 0)
		return;

	// Copy to host
	T* host_p;
	cudaMallocHost(&host_p, pitch*height*width);
	cudaMemcpy(host_p, dev_p, pitch*height*width, cudaMemcpyDeviceToHost);

	// Make indirect arrays
	T** host_pp = new T*[height*width];
	for (size_t y = 0; y < height*width; ++y)
		host_pp[y] = reinterpret_cast<T*>(reinterpret_cast<uint8_t*>(host_p) + y*pitch);
	T*** host_ppp = new T**[width];
	for (size_t x = 0; x < width; ++x)
		host_ppp[x] = &host_pp[x * height];

	// Callback
	callback(host_ppp);

	// Copy back
	cudaMemcpy(dev_p, host_p, pitch*height*width, cudaMemcpyHostToDevice);

	delete[] host_ppp;
	delete[] host_pp;
	cudaFreeHost(host_p);
}

}} // namespace nbl::cuda

#endif // __CUDA_MEM_SCOPE_H_
