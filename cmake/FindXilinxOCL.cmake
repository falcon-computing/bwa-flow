add_custom_target(XilinxOCL)

if(DEFINED ENV{XILINX_SDX})
  set(Xilinx_DIR "$ENV{XILINX_SDX}")
  set(XilinxOCL_INCLUDE_DIRS "${Xilinx_DIR}/runtime/include/1_2")
  set(XilinxOCL_LIBRARY_DIRS "${Xilinx_DIR}/runtime/lib/x86_64")
  set(XilinxOCL_LIBRARIES "xilinxopencl")
  set(XilinxOCL_FOUND TRUE)

  message(STATUS "Found Xilinx OpenCL Environment")
endif() 
