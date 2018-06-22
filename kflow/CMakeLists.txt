set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++0x -m64")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")


file(GLOB_RECURSE SRC_LIST src/*.cpp)
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/kflow
  ${Google_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS})

link_directories(
  ${Google_LIBRARY_DIRS}
  ${Boost_LIBRARY_DIRS})

add_library(kflow STATIC ${SRC_LIST})
add_dependencies(kflow Google)
add_custom_command(TARGET kflow POST_BUILD
  COMMAND "mv"             ARGS  libkflow.a ~libkflow.a
  COMMAND "${CMAKE_AR}"    ARGS  rs libkflow.a CMakeFiles/kflow.dir/src/*.o
  WORKING_DIRECTORY "${CMAKE_BINARY_DIRECTORY}")
