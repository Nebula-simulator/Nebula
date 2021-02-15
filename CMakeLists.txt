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

	if(CMAKE_CUDA_COMPILER_VERSION LESS 8.0)
		message(STATUS "Compiling without CUDA support: version 8.0 or higher required.")
		set(USE_CUDA OFF)
	else()
		set(USE_CUDA ON)

		set(CMAKE_CUDA_STANDARD 11)
		set(CMAKE_CUDA_STANDARD_REQUIRED ON)
		set(CMAKE_CUDA_EXTENSIONS OFF)
		add_definitions(-DCUDA_HEADERS_AVAILABLE=1)
		include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
	endif()
else()
	message(STATUS "Compiling without CUDA support: no compiler found")
	set(USE_CUDA OFF)
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
if(CMAKE_CUDA_COMPILER_VERSION LESS 11.0)
	string(APPEND CMAKE_CUDA_FLAGS " -gencode=arch=compute_30,code=[sm_30,compute_30]")
endif()
string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_35,code=[sm_35,compute_35]")
string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_37,code=[sm_37,compute_37]")
string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_50,code=[sm_50,compute_50]")
string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_52,code=[sm_52,compute_52]")
if(NOT CMAKE_CUDA_COMPILER_VERSION LESS 8.0)
	string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_60,code=[sm_60,compute_60]")
	string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_61,code=[sm_61,compute_61]")
endif()
if(NOT CMAKE_CUDA_COMPILER_VERSION LESS 9.0)
	string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_70,code=[sm_70,compute_70]")
endif()
if(NOT CMAKE_CUDA_COMPILER_VERSION LESS 10.0)
	string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_75,code=[sm_75,compute_75]")
endif()
if(NOT CMAKE_CUDA_COMPILER_VERSION LESS 11.0)
	string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_80,code=[sm_80,compute_80]")
	string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_86,code=[sm_86,compute_86]")
endif()

include_directories(3rdparty)
add_subdirectory(source)
