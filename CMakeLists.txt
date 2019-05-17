cmake_minimum_required(VERSION 3.1)

project(mesh_generator LANGUAGES CXX C)

option(BUILD_DOCUMENTATION "Build documentation" false)
option(MDADAMS_BUILD "Build for Mike" true)
option(SCAN_CONVERSION_REAL_TYPE_USE_RATIONAL
       "Use rational type as the scan conversion real type" false)
option(USE_TRIANGULATION_HIERARCHY "Use triangulation hierarchy" true)
set(DEBUG_LEVEL 1 CACHE string "Debug level")

set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/build/cmake/modules")

# This include should be placed as early as possible.
include(InSourceBuild)
include(Sanitizers)

set(CMAKE_VERBOSE_MAKEFILE true)

#include_directories(src/include)

if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
	  "Select the type of build (e.g., Debug, Release, RelWithDebInfo, MinSizeRel)."
	  FORCE)
endif()
message("CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED true)
set(CMAKE_VERBOSE_MAKEFILE true)

message("DEBUG_LEVEL ${DEBUG_LEVEL}")
message("ENABLE_ASSERTIONS ${ENABLE_ASSERTIONS}")
if (DEBUG_LEVEL GREATER 1)
	message(WARNING "Code will run slow with DEBUG_LEVEL > 1")
endif()
if (ENABLE_ASSERTIONS OR (DEBUG_LEVEL GREATER 0))
	include(EnableAssertions)
	message("Enabling assertions.")
endif()
if (DEBUG_LEVEL LESS 0)
	add_compile_options("-DNDEBUG")
endif()
add_definitions("-DDEBUG_LEVEL=${DEBUG_LEVEL}")

message("USE_TRIANGULATION_HIERARCHY ${USE_TRIANGULATION_HIERARCHY}")
if (USE_TRIANGULATION_HIERARCHY)
	add_definitions("-DUSE_TRIANGULATION_HIERARCHY")
endif()

message("SCAN_CONVERSION_REAL_TYPE_USE_RATIONAL ${SCAN_CONVERSION_REAL_TYPE_USE_RATIONAL}")
if (SCAN_CONVERSION_REAL_TYPE_USE_RATIONAL)
        add_definitions("-DSCAN_CONVERSION_REAL_TYPE_USE_RATIONAL")
endif()

if (MDADAMS_BUILD)
	set(SPL_TOP_DIR "/home/mdadams/local")
	set(SPLEL_DIR "${SPL_TOP_DIR}/lib64/cmake/SPLEL")
	set(SPL_DIR "${SPL_TOP_DIR}/lib64/cmake/SPL")
endif()

message("SPL_DIR ${SPL_DIR}")
message("SPLEL_DIR ${SPLEL_DIR}")

# Note: Performing find_package on CGAL appears to overwrite the value of
# Boost_LIBRARIES variable and possibly other Boost_* variables.
# This will cause problems if performing find_package on Boost is
# not performed after this.
find_package(CGAL REQUIRED)
message("Boost_LIBRARIES ${Boost_LIBRARIES}")
find_package(SPL REQUIRED)
find_package(SPLEL REQUIRED)
# Perform find_package on Boost after having performed find_package on CGAL
# (for reason stated above).
find_package(Boost COMPONENTS program_options REQUIRED)
message("Boost_LIBRARIES ${Boost_LIBRARIES}")

find_package(Doxygen)
find_package(LATEX)

find_program(SH_PROGRAM sh)
message("SH_PROGRAM ${SH_PROGRAM}")

add_compile_options("-g")

add_subdirectory(src)
#add_subdirectory(src/test)
add_subdirectory(src/helper)
if (BUILD_DOCUMENTATION)
	add_subdirectory(doc)
endif()

configure_file(bin/test_wrapper.in bin/test_wrapper)

if (SCAN_CONVERSION_REAL_TYPE_USE_RATIONAL)
	set(test_programs
  	  ${CMAKE_CURRENT_SOURCE_DIR}/src/test/run_test_rasterization
  	  ${CMAKE_CURRENT_SOURCE_DIR}/src/test/run_test_mesh_generator
	)
else ()
	set(test_programs
  	  ${CMAKE_CURRENT_SOURCE_DIR}/src/test/run_test_mesh_generator
	)
endif()
message("test_program ${test_programs}\n")
if (SH_PROGRAM)
        message("adding tests")
        foreach (program ${test_programs})
                get_filename_component(file "${program}" NAME)
                set(target "${file}")
                add_test(${target}
                  "${SH_PROGRAM}"
                  "${CMAKE_BINARY_DIR}/bin/test_wrapper"
                  "${program}"
                )
        endforeach()
endif()
enable_testing()


