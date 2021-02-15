add_subdirectory(io)
add_subdirectory(geometry/octree)

if(USE_CUDA)
	add_executable(nebula_gpu gpu_main.cu
		$<TARGET_OBJECTS:io>
		$<TARGET_OBJECTS:octree>)
	target_link_libraries(nebula_gpu ${LIBRARIES})
endif()

add_executable(nebula_cpu_mt cpu_mt_main.cpp
		$<TARGET_OBJECTS:io>
		$<TARGET_OBJECTS:octree>)
add_executable(nebula_cpu_edep cpu_energydep.cpp
		$<TARGET_OBJECTS:io>
		$<TARGET_OBJECTS:octree>)
target_link_libraries(nebula_cpu_mt ${LIBRARIES} Threads::Threads)
target_link_libraries(nebula_cpu_edep ${LIBRARIES} Threads::Threads)
