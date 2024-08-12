# This file is part of xhcfflib.
# SPDX-Identifier: LGPL-3.0-or-later
#
# xhcfflib is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# xhcfflib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

cmake_minimum_required(VERSION 3.14)

project(
  "xhcff"
  LANGUAGES "Fortran" "C" "CXX"
  VERSION "0.0.1"
  DESCRIPTION "Hydrostatic pressure force-field library"
)


# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# General configuration information
set(libs)
add_subdirectory("config")

###############################################################################
####################### SUBPROJECTS & DEPENDENCIES ############################
###############################################################################

# OpenMP parallelization
if(NOT TARGET "OpenMP::OpenMP_Fortran")
  if (WITH_OpenMP)
    find_package("OpenMP" REQUIRED)
    add_compile_definitions(WITH_OpenMP) 
  else()
    # Create dummy library
    add_library("OpenMP::OpenMP_Fortran" INTERFACE IMPORTED)
  endif()
endif()

# Linear algebra (contains BLAS)
if(NOT TARGET "LAPACK::LAPACK")
  find_package("LAPACK" REQUIRED)
endif()

# Fortran unit test interface (also used by other subprojects)
if(NOT TARGET "test-drive::test-drive" AND WITH_TESTS)
   find_package("test-drive" REQUIRED)
endif()

set(
  lib-deps
  "OpenMP::OpenMP_Fortran"
  "LAPACK::LAPACK"
  "$<$<VERSION_LESS:${CMAKE_VERSION},3.20>:BLAS::BLAS>"
)

# xhcff-lib sources
set(srcs)
add_subdirectory("src")


###############################################################################
############################### Static Library ################################
###############################################################################

# tight-binding library target
add_library(
  "${PROJECT_NAME}-lib"
  "${srcs}"
)
set_target_properties(
  "${PROJECT_NAME}-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  WINDOWS_EXPORT_ALL_SYMBOLS TRUE
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)
target_link_libraries(
  "${PROJECT_NAME}-lib"
  PUBLIC
  "${lib-deps}"
)
target_include_directories(
  "${PROJECT_NAME}-lib"
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)
if(NOT EXISTS "${PROJECT_BINARY_DIR}/include")
  file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/include")
endif()

###############################################################################
##################### Install and Export Instructions #########################
###############################################################################
# Export targets for other projects
add_library("${PROJECT_NAME}" INTERFACE)
target_link_libraries("${PROJECT_NAME}" INTERFACE "${PROJECT_NAME}-lib")
install(
  TARGETS
  "${PROJECT_NAME}"
  "${PROJECT_NAME}-lib"
  EXPORT
  "${PROJECT_NAME}-targets"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(
  EXPORT
  "${PROJECT_NAME}-targets"
  NAMESPACE
  "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
)
# Package license files
install(
  FILES
  "COPYING"
  "COPYING.LESSER"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)


# build a minimal app that can read xyz files and perform XHCFF singlepoints
if (build_exe)
add_subdirectory("app")
endif()

###############################################################################
################################## Testing ####################################
###############################################################################
# add the testsuite
if (WITH_TESTS)
  enable_testing()
  add_subdirectory("test")
  # Include C test program
  add_subdirectory("c_test")
endif()

