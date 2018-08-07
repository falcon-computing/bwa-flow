add_custom_target(IntelAltera)

set(IntelAltera_DIR "$ENV{XILINX_SDX}")
execute_process(COMMAND aocl compile-config OUTPUT_VARIABLE IntelAltera_INCLUDE_DIRS)

execute_process(COMMAND aocl ldflags OUTPUT_VARIABLE IntelAltera_Link_Dir)
string(REPLACE " " ";" IntelAltera_Link_Dir ${IntelAltera_Link_Dir})
foreach (subflag ${IntelAltera_Link_Dir})
  string(SUBSTRING ${subflag} 2 -1 iflag)
  list(APPEND IntelAltera_LIBRARY_DIRS ${iflag})
endforeach()

execute_process(COMMAND aocl ldlibs  OUTPUT_VARIABLE IntelAltera_Link_Libs)
string(REPLACE " " ";" IntelAltera_Link_Libs ${IntelAltera_Link_Libs})
list(FILTER IntelAltera_Link_Libs INCLUDE REGEX "^-l.+$" )
foreach (sublib ${IntelAltera_Link_Libs})
  string(SUBSTRING ${sublib} 2 -1 ilib)
  list(APPEND IntelAltera_LIBRARIES ${ilib})
endforeach()

set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-as-needed")