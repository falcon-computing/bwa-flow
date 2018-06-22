include(ExternalProject)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

cmake_minimum_required(VERSION 2.8.12)

project(fcs-bwa-flow)

# set up dependency
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package( Boost 1.53.0 REQUIRED
              COMPONENTS system thread iostreams filesystem regex )
find_package( ZLIB         REQUIRED )
find_package( Threads      REQUIRED )
find_package( Google       REQUIRED )
find_package( FalconLM     REQUIRED )


if ( USE_HTSLIB )
  find_package( HTS        REQUIRED )
endif()

if ( USE_FPGA )
  if ( DEFINED ENV{XILINX_SDX} )
    set ( FPGA_TYPE "xlnx" )
    find_package( Xilinx REQUIRED )
    set( FPGA_INCLUDE_DIRS "${Xilinx_INCLUDE_DIRS}" )
    set( FPGA_LIBRARY_DIRS "${Xilinx_LIBRARY_DIRS}")
    set( FPGA_LIBRARIES "${Xilinx_LIBRARIES}")
  elseif ( DEFINED ENV{ALTERAOCLSDKROOT} )
    set ( FPGA_TYPE "intel" )
    find_package( IntelAltera REQUIRED )
    set( FPGA_INCLUDE_DIRS "${IntelAltera_INCLUDE_DIRS}" )
    set( FPGA_LIBRARY_DIRS "${IntelAltera_LIBRARY_DIRS}")
    set( FPGA_LIBRARIES "${IntelAltera_LIBRARIES}")
    else()
    message( SEND_ERROR "No OpenCL library set up. Build without using FPGA." )
    set( USE_FPGA Off )
  endif()
endif()


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

# set up debug/release mode
if (NOT CMAKE_BUILD_TYPE)
  # message(STATUS "No build type selected, default to Debug")
  # set(CMAKE_BUILD_TYPE "Debug")
  message(STATUS "No build type selected, default to Debug")
  set(CMAKE_BUILD_TYPE "Debug")
endif()

add_subdirectory(kflow)
add_subdirectory(bwa)

# set up for compilation
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -fPIC")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -DUSELICENSE")

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

if ( USE_HTSLIB )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_HTSLIB" )
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
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${Google_INCLUDE_DIRS}
  ${FalconLM_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS})

include_directories( ./kflow/include )

link_directories(
  ${Google_LIBRARY_DIRS}
  ${FalconLM_LIBRARY_DIRS}
  ${Boost_LIBRARY_DIRS})

if ( USE_HTSLIB )
  include_directories(${HTS_INCLUDE_DIRS})
  link_directories(${HTS_LIBRARY_DIRS})
endif()

if ( USE_FPGA )
  include_directories(${FPGA_INCLUDE_DIRS})
  link_directories(${FPGA_LIBRARY_DIRS})
endif()


file(GLOB_RECURSE BWA_SRC_LIST src/*.cpp)
list(REMOVE_ITEM BWA_SRC_LIST "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp"
                          "${CMAKE_CURRENT_SOURCE_DIR}/src/mpi_main.cpp"
                          "${CMAKE_CURRENT_SOURCE_DIR}/src/MPIChannel.cpp"
                          "${CMAKE_CURRENT_SOURCE_DIR}/src/FPGAPipeline.cpp"
                          "${CMAKE_CURRENT_SOURCE_DIR}/src/SWTask.cpp"
                          "${CMAKE_CURRENT_SOURCE_DIR}/src/smithwaterman.cpp"
                          "${CMAKE_CURRENT_SOURCE_DIR}/src/XCLAgent.cpp"
                          "${CMAKE_CURRENT_SOURCE_DIR}/src/IntelAgent.cpp")

if ( USE_FPGA )
  list(APPEND BWA_SRC_LIST "${CMAKE_CURRENT_SOURCE_DIR}/src/FPGAPipeline.cpp"
                       "${CMAKE_CURRENT_SOURCE_DIR}/src/SWTask.cpp")
  if ( ${FPGA_TYPE} MATCHES "xlnx" )
    list(APPEND BWA_SRC_LIST "${CMAKE_CURRENT_SOURCE_DIR}/src/XCLAgent.cpp")
  elseif ( ${FPGA_TYPE} MATCHES "intel" )
    list(APPEND BWA_SRC_LIST "${CMAKE_CURRENT_SOURCE_DIR}/src/IntelAgent.cpp")
  endif()
endif()

#if ( USE_FPGA )
#  if ( ${FPGA_TYPE} MATCHES "xlnx" )
#    add_dependencies(.bwa-objs Xilinx)
#  elseif ( ${FPGA_TYPE} MATCHES "intel" )
#    add_dependencies(.bwa-objs IntelAltera)
#  endif()
#endif()

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
add_executable(bwa ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp ${BWA_SRC_LIST})
add_dependencies(bwa Google FalconLM)
if ( USE_HTSLIB )
  add_dependencies(bwa HTS)
endif()

set (CMAKE_STATIC_LINKER_FLAGS "--as-needed")
target_link_libraries(bwa
  cbwa
  kflow 
  ${Boost_LIBRARIES} 
  ${FalconLM_LIBRARIES}
  ${Google_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  ${ZLIB_LIBRARIES}
  ${CMAKE_DL_LIBS}
  ${HTS_LIBRARIES}
  ${FPGA_LIBRARIES})

#enable_testing()
#add_subdirectory(test)
