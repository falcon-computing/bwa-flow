set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-function -DUSE_MALLOC_WRAPPERS -DHAVE_PTHREAD")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Wall")

if ( USE_HTSLIB )
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_HTSLIB" )
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
if ( USE_HTSLIB )
  include_directories(${HTS_INCLUDE_DIRS})
  link_directories(${HTS_LIBRARY_DIRS})
endif()

file(GLOB_RECURSE SRC_LIST *.c)
list(REMOVE_ITEM SRC_LIST "${CMAKE_CURRENT_SOURCE_DIR}/main.c"
                          "${CMAKE_CURRENT_SOURCE_DIR}/example.c"
                          "${CMAKE_CURRENT_SOURCE_DIR}/gprof-helper.c")
add_library(cbwa STATIC ${SRC_LIST})
set_target_properties(cbwa PROPERTIES OUTPUT_NAME bwa)
if ( USE_HTSLIB )
  add_dependencies(cbwa HTS)
endif()
add_custom_command(TARGET cbwa POST_BUILD
  COMMAND "mv"             ARGS  libbwa.a ~libbwa.a
  COMMAND "${CMAKE_AR}"    ARGS  -csru libbwa.a CMakeFiles/cbwa.dir/*.o
  WORKING_DIRECTORY "${CMAKE_BINARY_DIRECTORY}")
