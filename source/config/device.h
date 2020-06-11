#ifndef __DEVICE_H_
#define __DEVICE_H_

/**
 * \file config/device.h
 * \brief Global definitions for GPU code
 */

/**
 * \def CUDA_COMPILER_AVAILABLE
 * \brief Set to 1 if the CUDA compiler is available, and to zero otherwise.
 *
 * This is not a setting! We just check if the CUDA compiler is being used.
 */
#ifdef __CUDACC__
	#define CUDA_COMPILER_AVAILABLE 1
	#define CUDA_HEADERS_AVAILABLE 1
#else
	#define CUDA_COMPILER_AVAILABLE 0
#endif

/**
 * \def CUDA_HEADERS_AVAILABLE
 * \brief Set to 1 if the CUDA headers are available, and to zero otherwise.
 *
 * This is actually set by CMake.
 *
 * The reason for having this preprocessor statement is the following. CMake
 * will compile .cu files with nvcc, and .cpp files with gcc. But when compiling
 * for a GPU device, we do want to use CUDA's float3 etc in .cpp code. So we
 * need to include the CUDA runtime in that situation.
 *
 * The present file provides a fallback definition for CUDA_HEADERS_AVAILABLE if
 * it is not already defined, based on whether nvcc is being used.
 */
#ifndef CUDA_HEADERS_AVAILABLE
	#define CUDA_HEADERS_AVAILABLE 0
#endif

/**
 * \def CUDA_COMPILING
 * \brief Detect if we are compiling for CUDA right now.
 *
 * This is not a setting! It's just a check if CUDA is the target right now.
 */
#ifdef __CUDA_ARCH__
	#define CUDA_COMPILING 1
#else
	#define CUDA_COMPILING 0
#endif


// Load CUDA headers
#if CUDA_HEADERS_AVAILABLE
	#include <cuda_runtime.h>
#endif


/**
 * \def CPU
 * \brief Designates that a given function must be run on a CPU, and never on a GPU.
 */

/**
 * \def GPU
 * \brief Designates that a given function must be run on a GPU, and never on a CPU.
 */

/**
 * \def PHYSICS
 * \brief Designates that a given function may either be run on a CPU or a GPU.
 * This is common for the physics code, hence the name.
 */
#if CUDA_COMPILER_AVAILABLE
	#define CPU __host__                // CPU-only code
	#define GPU __device__              // GPU-only code
	#define PHYSICS __host__ __device__ // Physics are GPU and CPU code
#else // CUDA_COMPILER_AVAILABLE
	#define CPU
	#define GPU
	#define PHYSICS
#endif // CUDA_COMPILER_AVAILABLE

#endif // __DEVICE_H_
