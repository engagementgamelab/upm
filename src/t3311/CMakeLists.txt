set (libname "t3311")
set (libdescription "Module for the Comet System T3311")
set (module_src ${libname}.cxx)
set (module_hpp ${libname}.hpp)

if (MODBUS_FOUND)
  set (reqlibname "libmodbus")
  include_directories(${MODBUS_INCLUDE_DIRS})
  upm_module_init()
  upm_target_link_libraries(${libname} ${MODBUS_LIBRARIES})
endif ()
