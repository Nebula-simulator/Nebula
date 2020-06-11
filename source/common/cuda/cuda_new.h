#ifndef __CUDA_NEW_H_
#define __CUDA_NEW_H_

/**
 * \file common/cuda/cuda_new.h
 * \brief Utility functions for allocating memory on the GPU
 */

namespace nbl { namespace cuda {

/**
 * \brief Allocate 1D array on the GPU.
 *
 * \tparam     T           Datatype to allocate memory for
 * \param[out] destination This function fills this with a pointer to the newly
 *                         allocated memory.
 * \param[in]  N           Number of elements to allocate
 */
template<typename T>
CPU void cuda_new(T** destination, size_t N)
{
	cudaError_t cudaStatus = cudaMalloc(destination, N * sizeof(T));
	if (cudaStatus != cudaSuccess)
	{
		throw std::bad_alloc();
	}
}

/**
 * \brief Allocate 2D array on the GPU.
 *
 * \tparam     T           Datatype to allocate memory for
 * \param[out] destination This function fills this with a pointer to the newly
 *                         allocated memory.
 * \param[out] pitch       This function fills this with the pitch of the newly
 *                         allocated memory.
 * \param[in]  width       Width of desired memory block
 * \param[in]  height      Height of desired memory block
 */
template<typename T>
CPU void cuda_new_2D(T** destination, size_t* pitch, size_t width, size_t height)
{
	cudaError_t cudaStatus = cudaMallocPitch(destination, pitch, height*sizeof(T), width);
	if (cudaStatus != cudaSuccess)
	{
		throw std::bad_alloc();
	}
}

/**
 * \brief Allocate 3D array on the GPU.
 *
 * \tparam     T           Datatype to allocate memory for
 * \param[out] destination This function fills this with a pointer to the newly
 *                         allocated memory.
 * \param[out] pitch       This function fills this with the pitch of the newly
 *                         allocated memory.
 * \param[in]  width       Width of desired memory block
 * \param[in]  height      Height of desired memory block
 * \param[in]  depth       Depth of desired memory block
 */
template<typename T>
CPU void cuda_new_3D(T** destination, size_t* pitch, size_t width, size_t height, size_t depth)
{
	cudaError_t cudaStatus = cudaMallocPitch(destination, pitch, depth*sizeof(T), height*width);
	if (cudaStatus != cudaSuccess)
	{
		throw std::bad_alloc();
	}
}

}} // namespace nbl::cuda

#endif // __CUDA_NEW_H_
