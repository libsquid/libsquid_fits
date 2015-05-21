# -------------------------- LICENSE -----------------------------------
#
# This file is part of the LibSQUID software library.
#
# LibSQUID is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# LibSQUID is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with LibSQUID.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2014 James Wren and Los Alamos National Laboratory
#

cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    "${CMAKE_SOURCE_DIR}/cmake/Modules")

# Use sensible default for undefined LIB_SUFFIX
if (NOT DEFINED LIB_SUFFIX)
  if(NOT (CMAKE_SIZEOF_VOID_P EQUAL 4))
    set(LIB_SUFFIX 64)
  endif()
endif()
# Force LIB_SUFFIX availability in ccmake
set(LIB_SUFFIX ${LIB_SUFFIX} CACHE STRING "Library suffix, appended to library path.")

# Parse squidfits.c SQUIDFITS_VERSION
file(STRINGS "${CMAKE_SOURCE_DIR}/squidfits.c" SQUIDFITS_VERSION_STRING REGEX
     "^[ \t]*#define[ \t]+SQUIDFITS_VERSION[ \t]+\"[0-9]+\\.[0-9]+\\.[0-9]+\"[ \t]*$")
if (SQUIDFITS_VERSION_STRING STREQUAL "")
  message(FATAL_ERROR "SQUIDFITS_VERSION not found in squidfits.c")
endif()
string(REGEX REPLACE ".*\"([0-9]+)\\.([0-9]+)\\.([0-9]+)\".*" "\\1.\\2.\\3"
       SQUIDFITS_VERSION_STRING "${SQUIDFITS_VERSION_STRING}")
# Parse squidfits.c SQUIDFITS_RELEASE
file(STRINGS "${CMAKE_SOURCE_DIR}/squidfits.c" SQUIDFITS_RELEASE_STRING REGEX
     "^[ \t]*#define[ \t]+SQUIDFITS_RELEASE[ \t]+\"[0-9]+\"[ \t]*$")
if (SQUIDFITS_RELEASE_STRING STREQUAL "")
  message(FATAL_ERROR "SQUIDFITS_RELEASE not found in squidfits.c")
endif()
string(REGEX REPLACE ".*\"([0-9]+)\".*" "\\1"
       SQUIDFITS_RELEASE_STRING "${SQUIDFITS_RELEASE_STRING}")

# Overall project options
project(SQUIDFITS)
set(PACKAGE_VERSION "${SQUIDFITS_VERSION_STRING}")
set(SQUIDFITS_VERSION "${SQUIDFITS_VERSION_STRING}")
set(SQUIDFITS_RELEASE "${SQUIDFITS_RELEASE_STRING}")

# Find necessary libraries
#libsquid_wcs
find_package(LIBSQUIDWCS REQUIRED)
if (${LIBSQUIDWCS_VERSION_STRING} VERSION_LESS "0.6.1")
  message(FATAL_ERROR "LIBSQUID_WCS version ${LIBSQUIDWCS_VERSION_STRING} is too old")
endif()
#libsquid
find_package(LIBSQUID REQUIRED)
if (${LIBSQUID_VERSION_STRING} VERSION_LESS "0.5.3")
  message(FATAL_ERROR "LIBSQUID version ${LIBSQUID_VERSION_STRING} is too old")
endif()
#cfitsio
find_package(CFITSIO REQUIRED)
if (${CFITSIO_VERSION_STRING} VERSION_LESS "3.0")
  message(FATAL_ERROR "CFITSIO version ${CFITSIO_VERSION_STRING} is too old")
endif()
#wcslib
find_package(WCSLIB REQUIRED)
if (${WCSLIB_VERSION_STRING} VERSION_LESS "4.19")
  message(FATAL_ERROR "WCSLIB version ${WCSLIB_VERSION_STRING} is too old")
endif()

# Build options
set(LIBS_PRIVATE
    ${LIBSQUIDWCS_LIBRARIES}
    ${LIBSQUID_LIBRARIES}
    ${CFITSIO_LIBRARIES}
    ${WCSLIB_LIBRARIES}
    m)
include_directories(
  ${LIBSQUIDWCS_INCLUDEDIR}
  ${LIBSQUID_INCLUDEDIR}
  ${CFITSIO_INCLUDEDIR}
  ${WCSLIB_INCLUDEDIR})
link_directories(${SQUIDFITS_BINARY_DIR})
link_libraries(${LIBS_PRIVATE})

# Generate list of source and header files
file(GLOB SQUIDFITS_HEADERS *.h)
file(GLOB SQUIDFITS_SOURCES *.c)

# Define targets and installs
foreach(INFILE ${SQUIDFITS_SOURCES})
  string(REGEX MATCH "[^/\\]+\\.c$" OUTFILE "${INFILE}")
  string(REPLACE ".c" "" OUTFILE ${OUTFILE})
  message("adding target ${OUTFILE} ${INFILE}")
  set(SQUIDFITS_TARGETS "${SQUIDFITS_TARGETS} ${OUTFILE}")
  add_executable(${OUTFILE} ${INFILE})
  install(TARGETS ${OUTFILE}
          DESTINATION bin
          COMPONENT squidfits)
endforeach(INFILE)

# Define packages
set(CPACK_PACKAGE_NAME "squidfits")
set(CPACK_PACKAGE_VERSION ${SQUIDFITS_VERSION})
set(CPACK_PACKAGE_RELEASE ${SQUIDFITS_RELEASE})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Programs to use LibSQUID library to manipulate FITS images")
set(CPACK_PACKAGE_CONTACT "Jim Wren, jwren@lanl.gov, Los Alamos National Laboratory")
set(CPACK_PACKAGE_VENDOR ${CPACK_PACKAGE_CONTACT})
set(CPACK_RPM_PACKAGE_REQUIRES
    "libsquid >= 0.5.3, libsquid_wcs >= 0.6.1, cfitsio >= 3.0, wcslib >= 4.19")
set(CPACK_DEBIAN_PACKAGE_DEPENDS
    "libsquid (>= 0.5.3), libsquid_wcs (>= 0.6.1) cfitsio (>= 3.0), wcslib (>= 4.19)")
set(CPACK_GENERATOR "RPM" "DEB")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "~$" "[.]bak$" "[.]swp$" "[.]svn$")
set(CPACK_PACKAGE_DIRECTORY ${CMAKE_BINARY_DIR}/packages)
set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(CPACK_RPM_EXCLUDE_FROM_AUTO_FILELIST
    "${CMAKE_INSTALL_PREFIX}"
    "${CMAKE_INSTALL_PREFIX}/bin")
set(CPACK_RPM_PACKAGE_RELEASE ${SQUIDFITS_RELEASE})
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_PACKAGE_RELEASE}.${CMAKE_SYSTEM_PROCESSOR}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_PACKAGE_RELEASE}")
include(CPack)

# Define package_source
add_custom_target(squidfits_dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
