ExternalProject_Add(glog-download
    PREFIX "deps/glog"
    URL https://s3.amazonaws.com/fcs-build-public/glog-falcon.tar.gz
    URL_MD5 2b1bb4285ef4c8963d5e0e338f1952b8
    SOURCE_DIR "${CMAKE_BINARY_DIR}/deps/glog/install"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

ExternalProject_Add(gflags-download
    PREFIX "deps/gflags"
    URL https://s3.amazonaws.com/fcs-build-public/gflags.tar.gz
    URL_MD5 1de8187489fbced5cc86c2ba241440e4
    SOURCE_DIR "${CMAKE_BINARY_DIR}/deps/gflags/install"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

ExternalProject_Add(googletest-download
    PREFIX "deps/googletest"
    URL https://s3.amazonaws.com/fcs-build-public/googletest.tar.gz
    URL_MD5 18fda945045354e264e3cca5428525d6
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
