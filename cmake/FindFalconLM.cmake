ExternalProject_Add(falconlm-download
    PREFIX "deps/falconlm"
    URL https://s3.amazonaws.com/fcs-build-public/falcon-lic.tar.gz
    URL_MD5 8f36e6cb08c7cb7a677d749bee1bfe6e
    CONFIGURE_COMMAND ""
    SOURCE_DIR "${CMAKE_BINARY_DIR}/deps/falconlm/install"
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

add_custom_target(FalconLM)
add_dependencies(FalconLM falconlm-download)

set(FalconLM_DIR "${CMAKE_BINARY_DIR}/deps/falconlm/install")
set(FalconLM_INCLUDE_DIRS "${FalconLM_DIR}/include")
set(FalconLM_LIBRARY_DIRS "${FalconLM_DIR}/lib")
set(FalconLM_LIBRARIES "falcon_license")
