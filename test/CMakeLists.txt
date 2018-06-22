if ( USE_HTSLIB )
  include_directories(${HTS_INCLUDE_DIRS})
  link_directories(${HTS_LIBRARY_DIRS})
endif()

if ( USE_FPGA )
  include_directories(${FPGA_INCLUDE_DIRS})
  link_directories(${FPGA_LIBRARY_DIRS})
endif()


file(GLOB_RECURSE BWA_TEST_LIST src/*.cpp)
list(REMOVE_ITEM BWA_TEST_LIST "${CMAKE_CURRENT_SOURCE_DIR}/src/FPGATests.cpp"
                          "${CMAKE_CURRENT_SOURCE_DIR}/src/ChannelTests.cpp")

if ( USE_FPGA )
  list(APPEND BWA_TEST_LIST "${CMAKE_CURRENT_SOURCE_DIR}/src/FPGATests.cpp")
endif()

add_executable(test_bwa ${BWA_TEST_LIST} ${BWA_SRC_LIST})

add_dependencies(test_bwa Google FalconLM)
if ( USE_HTSLIB )
  add_dependencies(test_bwa HTS)
endif()
#if ( USE_FPGA )
#  if ( ${FPGA_TYPE} MATCHES "xlnx" )
#    add_dependencies(test_bwa Xilinx)
#  elseif ( ${FPGA_TYPE} MATCHES "intel" )
#    add_dependencies(test_bwa IntelAltera)
#  endif()
#endif()

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

target_link_libraries(test_bwa
  kflow 
  cbwa
  ${Google_LIBRARIES}
  ${FPGA_LIBRARIES}
  ${HTS_LIBRARIES}
  ${FalconLM_LIBRARIES}
  ${Boost_LIBRARIES} 
  ${ZLIB_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  ${CMAKE_DL_LIBS})

add_test(bwa-test test_bwa)
