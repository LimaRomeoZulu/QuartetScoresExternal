# Genesis - A toolkit for working with phylogenetic data.
# Copyright (C) 2014-2017 Lucas Czech
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact:
# Lucas Czech <lucas.czech@h-its.org>
# Exelixis Lab, Heidelberg Institute for Theoretical Studies
# Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany

# --------------------------------------------------------------------------------------------------
#   CMake Init
# --------------------------------------------------------------------------------------------------

# Require the minimum version of any system that we test on - which is the version used by Travis.
# Thus, we can be sure that this works (as long as they don't change it...).
cmake_minimum_required (VERSION 3.0 FATAL_ERROR)
message (STATUS "CMake version ${CMAKE_VERSION}")

set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set (CMAKE_DISABLE_SOURCE_CHANGES  ON)

if (NOT WIN32)
    # Fun with colors!
    string(ASCII 27 Esc)

    set(ColorBlue   "${Esc}[34m") # Build information.
    set(ColorGreen  "${Esc}[32m") # Optional additions.
    set(ColorYellow "${Esc}[33m") # Tipps for the user.
    set(ColorRed    "${Esc}[31m") # Warnings and Errors.
    set(ColorEnd    "${Esc}[0m")
endif()

# --------------------------------------------------------------------------------------------------
#   Project Definition and Requirements
# --------------------------------------------------------------------------------------------------

project ( genesis CXX C)

# set (genesis_VERSION_MAJOR 0)
# set (genesis_VERSION_MINOR 1)

# Helpful output. Version information is unfortunately only available for CMake >= 2.8.9
message (STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} at ${CMAKE_CXX_COMPILER}")
message (STATUS "C compiler  : ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION} at ${CMAKE_C_COMPILER}")

# --------------------------------------------------------------------------------------------------
#   Build Options
# --------------------------------------------------------------------------------------------------

option (BUILD_SHARED_LIB    "Build shared library" ON)
option (BUILD_STATIC_LIB    "Build static library" OFF)

option (BUILD_APPLICATIONS  "Build applications"   OFF)
option (BUILD_PYTHON_MODULE "Build python module"  OFF)
option (BUILD_TESTS         "Build test suites"    OFF)

option (USE_THREADS         "Use a threading library (mainly, pthreads)"      OFF)
option (USE_OPENMP          "Use OpenMP if available"                         OFF)

# Define unity build option with list of valid values.
set( UNITY_BUILD FULL CACHE STRING "Build modules monolithically - trading memory for speed." )
set_property( CACHE UNITY_BUILD PROPERTY STRINGS FULL MODULES OFF )

if (BUILD_APPLICATIONS OR BUILD_PYTHON_MODULE OR BUILD_TESTS)
    # The applications, the python module and the tests all use the shared lib,
    # so we need to build it if one of those is activated.
    set (BUILD_SHARED_LIB ON)
endif()

# --------------------------------------------------------------------------------------------------
#   Build Type
# --------------------------------------------------------------------------------------------------

set (CMAKE_BUILD_TYPE DEBUG)

# Using the environment variable `GENESIS_DEBUG`, it is possible to change from
# release mode (default) to debug mode. This is useful for developers.
# It also activates test building automatically.
if( $ENV{GENESIS_DEBUG} )

    set (CMAKE_BUILD_TYPE DEBUG)

    # set (CMAKE_VERBOSE_MAKEFILE ON)
    SET( BUILD_TESTS    OFF )
    SET( UNITY_BUILD OFF )

endif()

message (STATUS "${ColorBlue}Building ${CMAKE_BUILD_TYPE}${ColorEnd}")
message (STATUS "Unity build: ${UNITY_BUILD}")

# --------------------------------------------------------------------------------------------------
#   Compiler and Linker Options
# --------------------------------------------------------------------------------------------------

# Set the warn flags to a very high level - except unknown pragmas, as this causes needless
# warnings with OpenMP and other pragma-based techniques.
set (WARN_FLAGS              "-Wall -Wextra -Wno-unknown-pragmas -pedantic -pedantic-errors")
set (CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -std=c++11 ${WARN_FLAGS}")

set (CMAKE_CXX_FLAGS_DEBUG   "-O2 -g -ggdb3 -DDEBUG -D_GLIBCXX_DEBUG")
set (CMAKE_CXX_FLAGS_RELEASE "-O2       -DNDEBUG")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Misc options
# STDCC = clang++
# MPICC   = mpic++
# mpi: CCFLAGS += -DUSE_MPI
# LDFLAGS = -lm $(THREADS)

# Pre-compilation check for endianness
include (TestBigEndian)
TEST_BIG_ENDIAN(IS_BIG_ENDIAN)
if( IS_BIG_ENDIAN )
    message(
        STATUS
        "${ColorYellow}Your system uses big endian memory, which is rare. Genesis is not tested "
        "with such systems. There might be bugs because of this!${ColorEnd}"
    )
endif()

# --------------------------------------------------------------------------------------------------
#   stxxl library
# --------------------------------------------------------------------------------------------------
# disallow in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# include the STXXL library
add_subdirectory(lib/stxxl/stxxl)
add_subdirectory(lib/genesis)

# apply STXXL CXXFLAGS to our configuration
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")

# add STXXL includes path
include_directories(${STXXL_INCLUDE_DIRS})
include_directories( ${PROJECT_SOURCE_DIR}/lib/genesis)

# build a program and link it with STXXL.
add_executable(QuartetScore QuartetScores.cpp)
target_link_libraries(QuartetScore ${STXXL_LIBRARIES})
target_link_libraries(QuartetScore genesis_lib_static)

set_property(TARGET QuartetScore PROPERTY CXX_STANDARD 11)
set_target_properties(QuartetScore PROPERTIES INSTALL_RPATH ${PROJECT_SOURCE_DIR}/build)

# --------------------------------------------------------------------------------------------------
#   Unity Build
# --------------------------------------------------------------------------------------------------

# Helper that builds monolithic source files by concatenating the original sources.
#
# It takes four parameters:
#
# input_directory: Where to look for source files (currently hardcoded to *.cpp)
# use_modules: Either "FULL" or "MODULES". The function offers two modes for monolithic building.
#     In "FULL" mode, all sources of the given input_directory are combined into one.
#     In "MODULES" mode, each immediate subdirectory of input_directory is treated separately,
#     resulting in one monolithic compilation unit per subdirectory. Files that are directly in
#     input_directory (not in its subdirectories) are simply used as normal single compilation units.
# name: Is used to distinguish between multiple monolithic parts during build (one for the lib,
#     one for test etc). The name needs to be uniq for each call of this function.
# resulting_source_files: Is used to store the resulting sources files. The value of this variable
#     can then ben used as source argument for add_library() and the like.
function( ASSEMBLE_MONOLITH input_directory use_modules name resulting_source_files  )

    # Helper: Return a list of all directories within a given directory, non-recursive.
    function( SUBDIRLIST result indir )
        file( GLOB children RELATIVE ${indir} ${indir}/* )
        set( dirlist "" )
        foreach(child ${children})
            if( IS_DIRECTORY ${indir}/${child} )
                LIST( APPEND dirlist ${child} )
            endif()
        endforeach()
        set( ${result} ${dirlist} PARENT_SCOPE )
    endfunction()

    # Helper: Append the contents of a file to another file.
    function( APPEND_TO_FILE IN_FILE OUT_FILE )
        file( READ ${IN_FILE} CONTENTS )
        file( APPEND ${OUT_FILE} "${CONTENTS}" )
    endfunction()

    # Clean previously assembled monolith files.
    file( REMOVE_RECURSE ${CMAKE_BINARY_DIR}/unity_assemble/${name} )

    IF( "${use_modules}" STREQUAL "MODULES" )

        # Process all modules (i.e., subdirs of the input dir).
        SUBDIRLIST( SUBDIRS ${input_directory} )
        foreach( SUBDIR ${SUBDIRS} )

            # Target files for this module.
            SET( ASSEMBLY_FILE ${CMAKE_BINARY_DIR}/unity_assemble/${name}/${SUBDIR}.cpp )
            SET( MONOLITH_FILE ${CMAKE_BINARY_DIR}/unity_sources/${name}/${SUBDIR}.cpp )

            # Find all sources of the module and concat them into a tmp assembly file.
            file (GLOB_RECURSE subdir_sources ${input_directory}/${SUBDIR}/*.cpp)
            file( WRITE ${ASSEMBLY_FILE} "" )
            foreach( subdir_source ${subdir_sources} )
                APPEND_TO_FILE( ${subdir_source} ${ASSEMBLY_FILE} )
            endforeach()

            # Copy the assemlby file to its destination. Using configure_file here ensures that we
            # only copy if the file changes, so that incremental builds are faster.
            configure_file( ${ASSEMBLY_FILE} ${MONOLITH_FILE} COPYONLY )
        endforeach()

        # Find all resulting source files. We need to treat top-level files in the lib dir
        # extra, as they are not a part of any module.
        file (GLOB top_level_sources ${input_directory}/*.cpp)
        file (GLOB unity_sources     ${CMAKE_BINARY_DIR}/unity_sources/${name}/*.cpp)
        set ( ${resulting_source_files} ${top_level_sources} ${unity_sources} PARENT_SCOPE )

    ELSEIF( "${use_modules}" STREQUAL "FULL" )

        # Target files.
        SET( ASSEMBLY_FILE ${CMAKE_BINARY_DIR}/unity_assemble/${name}/all.cpp )
        SET( MONOLITH_FILE ${CMAKE_BINARY_DIR}/unity_sources/${name}/all.cpp )

        # Find all sources and concat them into a tmp assembly file.
        file (GLOB_RECURSE all_sources ${input_directory}/*.cpp)
        file( WRITE ${ASSEMBLY_FILE} "" )
        foreach( subdir_source ${all_sources} )
            APPEND_TO_FILE( ${subdir_source} ${ASSEMBLY_FILE} )
        endforeach()

        # Copy the assemlby file to its destination. Using configure_file here ensures that we
        # only copy if the file changes, so that incremental builds are faster.
        configure_file( ${ASSEMBLY_FILE} ${MONOLITH_FILE} COPYONLY )

        # Return the assembled file.
        set ( ${resulting_source_files} ${MONOLITH_FILE} PARENT_SCOPE )

    ELSE()
        # Can only happen if some CMake file misuses this function.
        message( FATAL_ERROR
            "${ColorRed}ASSEMBLE_MONOLITH called with wrong arg: ${use_modules}${ColorEnd}"
        )
    ENDIF()
endfunction()

#----------------------------------------------------------------------------------------
