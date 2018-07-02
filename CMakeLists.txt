include(ExternalProject)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

cmake_minimum_required(VERSION 2.8.12)

project(bwa-flow)


# set up dependency
find_package( Boost 1.53.0 REQUIRED
              COMPONENTS system thread iostreams filesystem regex )
find_package( ZLIB         REQUIRED )
find_package( Threads      REQUIRED )

if ( NOT DEFINED USE_MPI )
  set( USE_MPI Off )
endif()
if ( USE_MPI )
  find_package( MPI        REQUIRED )
endif()

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package( Google       REQUIRED )
find_package( FalconLM     REQUIRED )
find_package( HTS        REQUIRED )

if ( DEFINED ENV{XILINX_SDX} )
  set( USE_FPGA On )
  set( FPGA_TYPE "xlnx" )
  find_package( Xilinx REQUIRED )
  set( FPGA_INCLUDE_DIRS "${Xilinx_INCLUDE_DIRS}" )
  set( FPGA_LIBRARY_DIRS "${Xilinx_LIBRARY_DIRS}")
  set( FPGA_LIBRARIES "${Xilinx_LIBRARIES}")
elseif ( DEFINED ENV{ALTERAOCLSDKROOT} )
  set( USE_FPGA On )
  set( FPGA_TYPE "intel" )
  find_package( IntelAltera REQUIRED )
  set( FPGA_INCLUDE_DIRS "${IntelAltera_INCLUDE_DIRS}" )
  set( FPGA_LIBRARY_DIRS "${IntelAltera_LIBRARY_DIRS}")
  set( FPGA_LIBRARIES "${IntelAltera_LIBRARIES}")
else()
  message( STATUS "No OpenCL module loaded. Build without using FPGA." )
  set( USE_FPGA Off )
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
  if ( ${FPGA_TYPE} MATCHES "xlnx" )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXILINX_FPGA" )
  elseif ( ${FPGA_TYPE} MATCHES "intel" )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DINTEL_FPGA" )
  endif()
endif()

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  set(BIN_VERSION "${BIN_VERSION}-dev")
elseif (${CMAKE_BUILD_TYPE} MATCHES "Release")
  if (DEFINED ENV{DEPLOYMENT})
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DDEPLOY_$ENV{DEPLOYMENT}")
    set(BIN_VERSION "${BIN_VERSION}-$ENV{DEPLOYMENT}")
  endif()
endif()

message(STATUS "version: ${BIN_VERSION}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERSION=\\\"${BIN_VERSION}\\\"")

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/src)

include_directories( ./kflow/include )

include_directories(
  ${Google_INCLUDE_DIRS}
  ${FalconLM_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${HTS_INCLUDE_DIRS})

link_directories(
  ${Google_LIBRARY_DIRS}
  ${FalconLM_LIBRARY_DIRS}
  ${Boost_LIBRARY_DIRS}
  ${HTS_LIBRARY_DIRS})

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


add_executable(bwa ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp ${BWA_SRC_LIST})
add_dependencies(bwa Google FalconLM HTS)

set (CMAKE_STATIC_LINKER_FLAGS "--as-needed")
target_link_libraries(bwa
  bwa-c
  kflow 
  ${Boost_LIBRARIES} 
  ${FalconLM_LIBRARIES}
  ${Google_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  ${ZLIB_LIBRARIES}
  ${CMAKE_DL_LIBS}
  ${HTS_LIBRARIES}
  ${FPGA_LIBRARIES})

if ( USE_MPI )
  add_executable(bwa-mpi ${CMAKE_CURRENT_SOURCE_DIR}/src/mpi_main.cpp ${BWA_MPI_SRC_LIST} ${BWA_SRC_LIST})
  add_dependencies(bwa-mpi Google FalconLM HTS)

  set (CMAKE_STATIC_LINKER_FLAGS "--as-needed")
  target_link_libraries(bwa-mpi
    bwa-c
    kflow 
    ${MPI_CXX_LIBRARIES}
    ${Boost_LIBRARIES} 
    ${FalconLM_LIBRARIES}
    ${Google_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${ZLIB_LIBRARIES}
    ${CMAKE_DL_LIBS}
    ${HTS_LIBRARIES}
    ${FPGA_LIBRARIES})
endif()

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   message(STATUS "Install prefix is set to ~/.falcon-genome/fcs-genome/${BIN_VERSION}")
   set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.falcon-genome/fcs-genome/${BIN_VERSION}")
endif()
install(TARGETS bwa RUNTIME DESTINATION .)
if ( USE_MPI )
  install(TARGETS bwa-mpi RUNTIME DESTINATION .)
endif()

enable_testing()
add_subdirectory(test)
