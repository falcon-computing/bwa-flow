ExternalProject_Add(falconlm-download
    PREFIX "falconlm"
    URL ${CMAKE_SOURCE_DIR}/deps/falcon-lic-latest.tgz
    CONFIGURE_COMMAND ""
    SOURCE_DIR "${CMAKE_BINARY_DIR}/falconlm/install"
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

add_custom_target(FalconLM)
add_dependencies(FalconLM falconlm-download)

set(FalconLM_DIR "${CMAKE_BINARY_DIR}/falconlm/install")
set(FalconLM_INCLUDE_DIRS "${FalconLM_DIR}/include")
set(FalconLM_LIBRARY_DIRS "${FalconLM_DIR}/lib")
set(FalconLM_LIBRARIES "falcon_license")
