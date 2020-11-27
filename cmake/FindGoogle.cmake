ExternalProject_Add(glog-download
    PREFIX "deps/glog"
    URL ${DEPS}/glog-falcon.tar.gz
    SOURCE_DIR "${CMAKE_BINARY_DIR}/deps/glog/install"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

ExternalProject_Add(gflags-download
    PREFIX "deps/gflags"
    URL ${DEPS}/gflags.tar.gz
    SOURCE_DIR "${CMAKE_BINARY_DIR}/deps/gflags/install"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

ExternalProject_Add(googletest-download
    PREFIX "deps/googletest"
    URL ${DEPS}/googletest.tar.gz
    SOURCE_DIR "${CMAKE_BINARY_DIR}/deps/googletest/install"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

add_custom_target(Google)
add_dependencies(Google
    glog-download
    gflags-download
    googletest-download)

set(Glog_DIR          ${CMAKE_BINARY_DIR}/deps/glog/install)
set(Glog_INCLUDE_DIRS ${Glog_DIR}/include)
set(Glog_LIBRARY_DIRS ${Glog_DIR}/lib)
set(Glog_LIBRARIES    glog)

set(Gflags_DIR          ${CMAKE_BINARY_DIR}/deps/gflags/install)
set(Gflags_INCLUDE_DIRS ${Gflags_DIR}/include)
set(Gflags_LIBRARY_DIRS ${Gflags_DIR}/lib)
set(Gflags_LIBRARIES    gflags)

set(Gtest_DIR          ${CMAKE_BINARY_DIR}/deps/googletest/install)
set(Gtest_INCLUDE_DIRS ${Gtest_DIR}/include)
set(Gtest_LIBRARY_DIRS ${Gtest_DIR}/lib)
set(Gtest_LIBRARIES    gtest)

set(Google_INCLUDE_DIRS
    ${Glog_INCLUDE_DIRS}
    ${Gflags_INCLUDE_DIRS}
    ${Gtest_INCLUDE_DIRS})

set(Google_LIBRARY_DIRS
    ${Glog_LIBRARY_DIRS}
    ${Gflags_LIBRARY_DIRS}
    ${Gtest_LIBRARY_DIRS})

set(Google_LIBRARIES
    ${Glog_LIBRARIES}
    ${Gflags_LIBRARIES}
    ${Gtest_LIBRARIES})
