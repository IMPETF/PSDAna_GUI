cmake_minimum_required(VERSION 2.8)

project(PSDAna_GUI)
set(PSDAna_GUI_VERSION_MAJOR 1)
set(PSDAna_GUI_VERSION_MINOR 1)

#system configuration
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PSDAna_GUI_SOURCE_DIR}/cmake)
#message(STATUS "${CMAKE_MODULE_PATH}")
#include(CMakeParseArguments)
#set ROOT system
if(NOT DEFINED ENV{ROOTSYS})
    message(FATAL_ERROR "Please set environment variable ROOTSYS in your system")
else()
    set(ROOTSYS $ENV{ROOTSYS})
    message(STATUS "ROOTSYS = ${ROOTSYS}")
endif()

#root executbles
if(WIN32)
    find_program(ROOTCINT_EXECUTABLE rootcint PATHS ${ROOTSYS}/bin)
elseif(UNIX)
    find_program(ROOTCINT_EXECUTABLE rootcint PATHS ${ROOTSYS}/bin)
    find_program(ROOT_CONFIG_EXECUTABLE root-config PATHS ${ROOTSYS}/bin)
else(WIN32)
    message(FATAL_ERROR "Platform not supportted: only windows and unix are supported")
endif(WIN32)

#root include path
#set(ROOT_INCLUDE_DIR)
if(WIN32)
    set(ROOT_INCLUDE_DIR ${ROOTSYS}/include)
else(WIN32)
    execute_process(
        COMMAND ${ROOT_CONFIG_EXECUTABLE} --incdir
        OUTPUT_VARIABLE ROOT_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endif(WIN32)
include_directories(${ROOT_INCLUDE_DIR})
message(STATUS "ROOT_INCLUDE_DIR = ${ROOT_INCLUDE_DIR}")

#root library
if(WIN32)
    set(ROOT_LIBRARY_DIR ${ROOTSYS}/lib)
    set(ROOT_LIBRARIES ${ROOTSYS}/lib/libCore.lib
                 ${ROOTSYS}/lib/libCint.lib ${ROOTSYS}/lib/libHist.lib
                 ${ROOTSYS}/lib/libGraf.lib ${ROOTSYS}/lib/libGraf3d.lib
                 ${ROOTSYS}/lib/libGpad.lib ${ROOTSYS}/lib/libTree.lib
                 ${ROOTSYS}/lib/libRint.lib ${ROOTSYS}/lib/libPostscript.lib
                 ${ROOTSYS}/lib/libMatrix.lib ${ROOTSYS}/lib/libPhysics.lib
                 ${ROOTSYS}/lib/libNet.lib ${ROOTSYS}/lib/libRIO.lib
                 ${ROOTSYS}/lib/libMathCore.lib
                 ${ROOTSYS}/lib/libSpectrum.lib
                 ${ROOTSYS}/lib/libGui.lib ${ROOTSYS}/lib/libGraf.lib
                 ${ROOTSYS}/lib/libGpad.lib)
else(WIN32)
    execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --libdir
        OUTPUT_VARIABLE ROOT_LIBRARY_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --libs
        OUTPUT_VARIABLE ROOT_LIBRARIES
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endif(WIN32)
link_directories(${ROOT_LIBRARY_DIR})
message(STATUS "ROOT_LIBRARY_DIR = ${ROOT_LIBRARY_DIR}")
message(STATUS "ROOT_LIBRARIES = ${ROOT_LIBRARIES}")

#process root dictionary
#define a function to process dictionary generation
include(CMakeParseArguments)
function(ROOT_GENERATE_DICTIONARY dictionary)
  CMAKE_PARSE_ARGUMENTS(ARG "" "" "LINKDEF;OPTIONS" "" ${ARGN})
  #---Get the list of header files-------------------------
  set(headerfiles)
  foreach(fp ${ARG_UNPARSED_ARGUMENTS})
    file(GLOB files ${fp})
    if(files)
      foreach(f ${files})
        if(NOT f MATCHES LinkDef)
          set(headerfiles ${headerfiles} ${f})
        endif()
      endforeach()
    else()
      set(headerfiles ${headerfiles} ${fp})
    endif()
  endforeach()
  #---Get the list of include directories------------------
  get_directory_property(incdirs INCLUDE_DIRECTORIES)
  set(includedirs)
  foreach( d ${incdirs})
   set(includedirs ${includedirs} -I${d})
  endforeach()
  #---Get LinkDef.h file------------------------------------
  set(linkdefs)
  foreach( f ${ARG_LINKDEF})
    if( IS_ABSOLUTE ${f})
      set(linkdefs ${linkdefs} ${f})
    else()
      if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/inc/${f})
        set(linkdefs ${linkdefs} ${CMAKE_CURRENT_SOURCE_DIR}/inc/${f})
      else()
        set(linkdefs ${linkdefs} ${CMAKE_CURRENT_SOURCE_DIR}/${f})
      endif()
    endif()
  endforeach()
  #---call rootcint------------------------------------------
  add_custom_command(OUTPUT ${dictionary}.cxx ${dictionary}.h
                     COMMAND ${ROOTCINT_EXECUTABLE} -f  ${dictionary}.cxx
                                          -c ${ARG_OPTIONS} ${includedirs} ${headerfiles} ${linkdefs}
                     DEPENDS ${headerfiles} ${linkdefs})
endfunction()
#add header files that need to generate dictionary
set(ROOT_DICT_SUFFIX Dict)
set(ROOT_LINKDEF_SUFFIX LinkDef.h)
set(ROOT_DICT_PRESOURCES GuiFrame)
set(ROOT_DICT_SOURCES)
foreach(name ${ROOT_DICT_PRESOURCES})
    ROOT_GENERATE_DICTIONARY(${name}${ROOT_DICT_SUFFIX} ${name}.h LINKDEF ${name}${ROOT_LINKDEF_SUFFIX})
    set(ROOT_DICT_SOURCES ${ROOT_DICT_SOURCES} ${name}${ROOT_DICT_SUFFIX}.cxx)
endforeach()


#add executable
set(PRIMARY_SOURCES main.cxx GuiFrame.cxx)
set(TOTAL_SOURCES ${PRIMARY_SOURCES} ${ROOT_DICT_SOURCES})
add_executable(PSDAna_GUI ${TOTAL_SOURCES})
target_link_libraries(PSDAna_GUI ${ROOT_LIBRARIES})

#install the targets
install(TARGETS PSDAna_GUI DESTINATION bin)
install(DIRECTORY documents/ DESTINATION doc)

#build a installer package
include(InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set (CPACK_PACKAGE_VERSION_MAJOR "${PSDAna_GUI_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${PSDAna_GUI_VERSION_MINOR}")
set (CPACK_PACKAGE_EXECUTABLES PSDAna_GUI "PSDAna_GUI")
include (CPack)
