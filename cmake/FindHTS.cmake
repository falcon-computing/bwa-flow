ExternalProject_Add(hts-download
    PREFIX "deps/hts"
    URL ${DEPS}/htslib-1.3.1.tar.gz
    CONFIGURE_COMMAND ""
    SOURCE_DIR "${CMAKE_BINARY_DIR}/deps/hts/install"
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

add_custom_target(HTS)
add_dependencies(HTS hts-download)

set(HTS_DIR "${CMAKE_BINARY_DIR}/deps/hts/install")
set(HTS_INCLUDE_DIRS "${HTS_DIR}/include")
set(HTS_LIBRARY_DIRS "${HTS_DIR}/lib")
set(HTS_LIBRARIES "hts")
