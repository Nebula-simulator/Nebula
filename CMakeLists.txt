cmake_minimum_required(VERSION 3.8)
project(nebula C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)


# Find CUDA, if available
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
	enable_language(CUDA)
	set(CMAKE_CUDA_STANDARD 11)
	set(CMAKE_CUDA_STANDARD_REQUIRED ON)
	set(CMAKE_CUDA_EXTENSIONS OFF)
	add_definitions(-DCUDA_HEADERS_AVAILABLE=1)
	include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
else()
	message(STATUS "Compiling without CUDA support")
endif()


# Find HDF5 library
find_package(HDF5 1.8.13 REQUIRED C HL)
include_directories(${HDF5_INCLUDE_DIRS})
add_definitions(${HDF5_DEFINITIONS})
list(APPEND LIBRARIES ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

# Find MPI library if compiling against parallel HDF5
if(HDF5_IS_PARALLEL)
	find_package(MPI REQUIRED)
	include_directories(${MPI_CXX_INCLUDE_PATH})
	list(APPEND LIBRARIES ${MPI_CXX_LIBRARIES})
endif()

# Find threading library
find_package(Threads REQUIRED)


if(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE "Release" CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
		FORCE)
endif()

string(APPEND CMAKE_CUDA_FLAGS " --expt-relaxed-constexpr")
string(APPEND CMAKE_CUDA_FLAGS " -prec-sqrt=false")
string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_30,code=compute_30")
string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_37,code=compute_37")
string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_60,code=compute_60")
string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_61,code=compute_61")
if(CMAKE_CUDA_COMPILER_VERSION LESS 9.0)
	string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_20,code=compute_20")
endif(CMAKE_CUDA_COMPILER_VERSION LESS 9.0)

include_directories(3rdparty)
add_subdirectory(source)
