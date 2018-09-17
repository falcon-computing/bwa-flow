ExternalProject_Add(falconlm-download
    PREFIX "falconlm"
    URL https://s3.amazonaws.com/fcs-build-public/falcon-lic-v1.2.tgz
    URL_MD5 911e5526bbbb4b2967580173b59cc0ff 
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
