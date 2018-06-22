ExternalProject_Add(hts-download
    PREFIX "deps/hts"
    URL https://s3.amazonaws.com/fcs-build-public/htslib-1.3.1.tar.gz
    URL_MD5 7d0273d4a447b8c01d1a82292eea1ad0
    CONFIGURE_COMMAND ""
    SOURCE_DIR "${CMAKE_BINARY_DIR}/deps/hts/install"
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

add_custom_target(HTS)
add_dependencies(HTS hts-download)

set(HTS_DIR "${CMAKE_BINARY_DIR}/deps/hts/install")
set(HTS_INCLUDE_DIRS "${HTS_DIR}")
set(HTS_LIBRARY_DIRS "${HTS_DIR}")
set(HTS_LIBRARIES "hts")
