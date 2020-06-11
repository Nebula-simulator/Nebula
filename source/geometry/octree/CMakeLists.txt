list(APPEND SOURCE_FILES octree_builder.cpp)
list(APPEND SOURCE_FILES tribox.cc)
add_library(octree OBJECT ${SOURCE_FILES})
