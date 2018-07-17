set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_MALLOC_WRAPPERS -DHAVE_PTHREAD -DUSE_HTSLIB")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wno-unused-function -O2")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Wall")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    ${HTS_INCLUDE_DIRS})
link_directories(${HTS_LIBRARY_DIRS})

file(GLOB_RECURSE SRC_LIST *.c)
list(REMOVE_ITEM SRC_LIST "${CMAKE_CURRENT_SOURCE_DIR}/main.c"
                          "${CMAKE_CURRENT_SOURCE_DIR}/example.c"
                          "${CMAKE_CURRENT_SOURCE_DIR}/gprof-helper.c")
add_library(bwa-c STATIC ${SRC_LIST})
set_target_properties(bwa-c PROPERTIES OUTPUT_NAME bwa)
add_dependencies(bwa-c HTS)
