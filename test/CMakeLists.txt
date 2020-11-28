file(GLOB_RECURSE BWA_TEST_LIST src/*.cpp)
list(REMOVE_ITEM BWA_TEST_LIST "${CMAKE_CURRENT_SOURCE_DIR}/src/FPGATests.cpp"
                          "${CMAKE_CURRENT_SOURCE_DIR}/src/ChannelTests.cpp")

if ( USE_FPGA )
  list(APPEND BWA_TEST_LIST "${CMAKE_CURRENT_SOURCE_DIR}/src/FPGATests.cpp")
endif()
if ( USE_MPI )
  list(APPEND BWA_TEST_LIST "${CMAKE_CURRENT_SOURCE_DIR}/src/ChannelTests.cpp")
  list(APPEND BWA_SRC_LIST "${BWA_MPI_SRC_LIST}")
endif()

add_executable(test_bwa ${BWA_TEST_LIST} ${BWA_SRC_LIST})

add_dependencies(test_bwa Google HTS)

target_link_libraries(test_bwa
  kflow 
  bwa-c
  ${MPI_CXX_LIBRARIES}
  ${FPGA_LIBRARIES}
  ${HTS_LIBRARIES}
  ${Google_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  ${ZLIB_LIBRARIES}
  ${CMAKE_DL_LIBS}
  ${Boost_LIBRARIES})

add_test(NAME    bwa-test
         COMMAND bash test-driver.sh $<TARGET_FILE:test_bwa>
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test)
