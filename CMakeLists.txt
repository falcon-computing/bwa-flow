include(ExternalProject)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

cmake_minimum_required(VERSION 2.8.12)

project(bwa-flow)

# set up dependency
find_package(Boost 1.53.0 REQUIRED
             COMPONENTS system thread chrono iostreams filesystem regex )
find_package(ZLIB         REQUIRED)
find_package(Threads      REQUIRED)

if ( NOT DEFINED USE_MPI )
  set( USE_MPI Off )
endif()
if ( USE_MPI )
  find_package( MPI        REQUIRED )
endif()

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(Google REQUIRED)
find_package(HTS REQUIRED)
find_package(XilinxOCL)
find_package(IntelOCL)

if(XilinxOCL_FOUND)
  set(USE_FPGA On )
  message(STATUS "Build with Xilinx FPGA")
  set(FPGA_TYPE "xlnx" )
  set(FPGA_INCLUDE_DIRS "${XilinxOCL_INCLUDE_DIRS}")
  set(FPGA_LIBRARY_DIRS "${XilinxOCL_LIBRARY_DIRS}")
  set(FPGA_LIBRARIES "${XilinxOCL_LIBRARIES}")

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXILINX_FPGA" )
elseif(IntelOCL_FOUND)
  set(USE_FPGA On )
  message(STATUS "Build with Intel FPGA")
  set(FPGA_TYPE "intel")
  set(FPGA_INCLUDE_DIRS "${IntelOCL_INCLUDE_DIRS}")
  set(FPGA_LIBRARY_DIRS "${IntelOCL_LIBRARY_DIRS}")
  set(FPGA_LIBRARIES "${IntelOCL_LIBRARIES}")

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DINTEL_FPGA" )
else()
  message(STATUS "No OpenCL module loaded. Build without using FPGA." )
  set(USE_FPGA Off)
endif()

# set up debug/release mode
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Debug")
  set(CMAKE_BUILD_TYPE "Debug")
endif()
message( STATUS "Build type: ${CMAKE_BUILD_TYPE}" )

add_subdirectory(kflow)
add_subdirectory(bwa bwa_c)

# set up for compilation
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -fPIC -DUSE_HTSLIB")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -DUSELICENSE")
if ( USE_MPI )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_LINK_FLAGS}")
endif()


# get code version
execute_process( COMMAND
  git describe --tags --always
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE )

set(BIN_VERSION "${GIT_VERSION}")
if ( USE_FPGA )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBUILD_FPGA" )
  set(BIN_VERSION "${BIN_VERSION}-${FPGA_TYPE}")
endif()

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  set(BIN_VERSION "${BIN_VERSION}-dev")
endif()
if (DEFINED DEPLOYMENT_DST)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEPLOY_${DEPLOYMENT_DST}")
  set(BIN_VERSION "${BIN_VERSION}-${DEPLOYMENT_DST}")
endif()

message(STATUS "version: ${BIN_VERSION}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERSION=\\\"${BIN_VERSION}\\\"")

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ./kflow/include
  ${Google_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${HTS_INCLUDE_DIRS})

link_directories(
  ${Google_LIBRARY_DIRS}
  ${Boost_LIBRARY_DIRS}
  ${HTS_LIBRARY_DIRS}
  ${ZLIB_LIBRARY_DIRS})

file(GLOB BWA_SRC_LIST src/*.cpp)
list(REMOVE_ITEM BWA_SRC_LIST "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
list(REMOVE_ITEM BWA_SRC_LIST "${CMAKE_CURRENT_SOURCE_DIR}/src/mpi_main.cpp")

if ( USE_FPGA )
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/fpga
                      ${CMAKE_CURRENT_SOURCE_DIR}/src/fpga/${FPGA_TYPE}
                      ${FPGA_INCLUDE_DIRS})
  link_directories(${FPGA_LIBRARY_DIRS})

  file(GLOB BWA_FPGA_SRC_LIST src/fpga/*.cpp)
  list(APPEND BWA_SRC_LIST ${BWA_FPGA_SRC_LIST})

  file(GLOB BWA_FPGA_EXTRA_SRC_LIST src/fpga/${FPGA_TYPE}/*.cpp)
  list(APPEND BWA_SRC_LIST ${BWA_FPGA_EXTRA_SRC_LIST})
endif()

if ( USE_MPI )
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/mpi
                      ${MPI_CXX_INCLUDE_PATH})
  file(GLOB BWA_MPI_SRC_LIST  src/mpi/*.cpp)
endif()

add_executable(bwa-flow 
  ${BWA_SRC_LIST}
  ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)

add_dependencies(bwa-flow Google HTS)

set (CMAKE_STATIC_LINKER_FLAGS "--as-needed")
target_link_libraries(bwa-flow
  bwa-c
  kflow 
  ${HTS_LIBRARIES}
  ${Google_LIBRARIES}
  ${Boost_LIBRARIES}
  ${ZLIB_LIBRARIES}
  ${CMAKE_DL_LIBS}
  ${FPGA_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT})

if ( USE_MPI )
  add_executable(bwa-mpi ${CMAKE_CURRENT_SOURCE_DIR}/src/mpi_main.cpp ${BWA_MPI_SRC_LIST} ${BWA_SRC_LIST})
  add_dependencies(bwa-mpi Google HTS)

  set (CMAKE_STATIC_LINKER_FLAGS "--as-needed")
  target_link_libraries(bwa-mpi
    bwa-c
    kflow 
    ${MPI_CXX_LIBRARIES}
    ${FPGA_LIBRARIES}
    ${HTS_LIBRARIES}
    ${Google_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${ZLIB_LIBRARIES}
    ${CMAKE_DL_LIBS}
    ${Boost_LIBRARIES})
endif()

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   message(STATUS "Install prefix is set to ~/.falcon-genome/bwa/${BIN_VERSION}")
   set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.falcon-genome/bwa/${BIN_VERSION}")
endif()
install(TARGETS bwa-flow RUNTIME DESTINATION .)
if ( USE_MPI )
  install(TARGETS bwa-mpi RUNTIME DESTINATION .)
endif()

enable_testing()
add_subdirectory(test)
