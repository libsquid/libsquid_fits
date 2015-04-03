# - Try to find LIBSQUIDWCS
# Once done this will define
#
#  LIBSQUIDWCS_FOUND - system has LIBSQUIDWCS
#  LIBSQUIDWCS_INCLUDE_DIR - the LIBSQUIDWCS include directory
#  LIBSQUIDWCS_LIBRARIES - Link these to use LIBSQUIDWCS
#  LIBSQUIDWCS_VERSION_STRING - Human readable version number of cfitsio
#  LIBSQUIDWCS_VERSION_MAJOR  - Major version number of cfitsio
#  LIBSQUIDWCS_VERSION_MINOR  - Minor version number of cfitsio
#
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

if (LIBSQUIDWCS_INCLUDE_DIR AND LIBSQUIDWCS_LIBRARIES)

  # in cache already, be quiet
  set(LIBSQUIDWCS_FIND_QUIETLY TRUE)

else (LIBSQUIDWCS_INCLUDE_DIR AND LIBSQUIDWCS_LIBRARIES)

  find_path(LIBSQUIDWCS_INCLUDE_DIR libsquid_wcs.h
    PATHS
    $ENV{LIBSQUIDWCS}
    ${CMAKE_SOURCE_DIR}/libsquid_wcs
    ${CMAKE_SOURCE_DIR}/../libsquid_wcs
    ${_obIncDir}
    ${GNUWIN32_DIR}/include
  )
  message("LIBSQUIDWCS_INCLUDE_DIR ${LIBSQUIDWCS_INCLUDE_DIR}")

  find_library(LIBSQUIDWCS_LIBRARIES NAMES squid_wcs
    PATHS
    $ENV{LIBSQUIDWCS}
    ${CMAKE_BINARY_DIR}/libsquid_wcs
    ${_obLinkDir}
    ${GNUWIN32_DIR}/lib
  )
  message("LIBSQUIDWCS_LIBRARIES ${LIBSQUIDWCS_LIBRARIES}")

  if(LIBSQUIDWCS_INCLUDE_DIR AND LIBSQUIDWCS_LIBRARIES)
    set(LIBSQUIDWCS_FOUND TRUE)
  else (LIBSQUIDWCS_INCLUDE_DIR AND LIBSQUIDWCS_LIBRARIES)
    set(LIBSQUIDWCS_FOUND FALSE)
  endif(LIBSQUIDWCS_INCLUDE_DIR AND LIBSQUIDWCS_LIBRARIES)

  if (LIBSQUIDWCS_FOUND)

    # Find the version of the header
    FILE(READ "${LIBSQUIDWCS_INCLUDE_DIR}/libsquid_wcs.h" LIBSQUIDWCS_H)
    STRING(REGEX REPLACE ".*#define LIBSQUIDWCS_VERSION[^0-9]*([0-9]+)\\.([0-9]+)\\.([0-9]+).*" "\\1.\\2.\\3" LIBSQUIDWCS_VERSION_STRING "${LIBSQUIDWCS_H}")
    STRING(REGEX REPLACE "^([0-9]+)[.]([0-9]+)[.]([0-9]+)" "\\1" LIBSQUIDWCS_VERSION_MAJOR ${LIBSQUIDWCS_VERSION_STRING})
    STRING(REGEX REPLACE "^([0-9]+)[.]([0-9]+)[.]([0-9]+)" "\\2" LIBSQUIDWCS_VERSION_MINOR ${LIBSQUIDWCS_VERSION_STRING})
    STRING(REGEX REPLACE "^([0-9]+)[.]([0-9]+)[.]([0-9]+)" "\\3" LIBSQUIDWCS_VERSION_PATCH ${LIBSQUIDWCS_VERSION_STRING})
    message(STATUS "found version string ${LIBSQUIDWCS_VERSION_STRING}")

    if (NOT LIBSQUIDWCS_FIND_QUIETLY)
      message(STATUS "Found LIBSQUIDWCS ${LIBSQUIDWCS_VERSION_MAJOR}.${LIBSQUIDWCS_VERSION_MINOR}.${LIBSQUIDWCS_VERSION_PATCH}: ${LIBSQUIDWCS_LIBRARIES}")
    endif (NOT LIBSQUIDWCS_FIND_QUIETLY)
  else (LIBSQUIDWCS_FOUND)
    if (LIBSQUIDWCS_FIND_REQUIRED)
      message(STATUS "LIBSQUIDWCS not found.")
    endif (LIBSQUIDWCS_FIND_REQUIRED)
  endif (LIBSQUIDWCS_FOUND)

  mark_as_advanced(LIBSQUIDWCS_INCLUDE_DIR LIBSQUIDWCS_LIBRARIES)

endif (LIBSQUIDWCS_INCLUDE_DIR AND LIBSQUIDWCS_LIBRARIES)
