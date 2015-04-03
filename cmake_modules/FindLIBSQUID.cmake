# - Try to find LIBSQUID
# Once done this will define
#
#  LIBSQUID_FOUND - system has LIBSQUID
#  LIBSQUID_INCLUDE_DIR - the LIBSQUID include directory
#  LIBSQUID_LIBRARIES - Link these to use LIBSQUID
#  LIBSQUID_VERSION_STRING - Human readable version number of cfitsio
#  LIBSQUID_VERSION_MAJOR  - Major version number of cfitsio
#  LIBSQUID_VERSION_MINOR  - Minor version number of cfitsio
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


if (LIBSQUID_INCLUDE_DIR AND LIBSQUID_LIBRARIES)

  # in cache already, be quiet
  set(LIBSQUID_FIND_QUIETLY TRUE)

else (LIBSQUID_INCLUDE_DIR AND LIBSQUID_LIBRARIES)

  find_path(LIBSQUID_INCLUDE_DIR libsquid.h
    PATHS
    $ENV{LIBSQUID}
    ${CMAKE_SOURCE_DIR}/libsquid
    ${CMAKE_SOURCE_DIR}/../libsquid
    ${_obIncDir}
    ${GNUWIN32_DIR}/include
  )
  message("LIBSQUID_INCLUDE_DIR ${LIBSQUID_INCLUDE_DIR}")

  find_library(LIBSQUID_LIBRARIES NAMES squid
    PATHS
    $ENV{LIBSQUID}
    ${CMAKE_BINARY_DIR}/libsquid
    ${_obLinkDir}
    ${GNUWIN32_DIR}/lib
  )
  message("LIBSQUID_LIBRARIES ${LIBSQUID_LIBRARIES}")

  if(LIBSQUID_INCLUDE_DIR AND LIBSQUID_LIBRARIES)
    set(LIBSQUID_FOUND TRUE)
  else (LIBSQUID_INCLUDE_DIR AND LIBSQUID_LIBRARIES)
    set(LIBSQUID_FOUND FALSE)
  endif(LIBSQUID_INCLUDE_DIR AND LIBSQUID_LIBRARIES)

  if (LIBSQUID_FOUND)

    # Find the version of the header
    FILE(READ "${LIBSQUID_INCLUDE_DIR}/libsquid.h" LIBSQUID_H)
    STRING(REGEX REPLACE ".*#define LIBSQUID_VERSION[^0-9]*([0-9]+)\\.([0-9]+)\\.([0-9]+).*" "\\1.\\2.\\3" LIBSQUID_VERSION_STRING "${LIBSQUID_H}")
    STRING(REGEX REPLACE "^([0-9]+)[.]([0-9]+)[.]([0-9]+)" "\\1" LIBSQUID_VERSION_MAJOR ${LIBSQUID_VERSION_STRING})
    STRING(REGEX REPLACE "^([0-9]+)[.]([0-9]+)[.]([0-9]+)" "\\2" LIBSQUID_VERSION_MINOR ${LIBSQUID_VERSION_STRING})
    STRING(REGEX REPLACE "^([0-9]+)[.]([0-9]+)[.]([0-9]+)" "\\3" LIBSQUID_VERSION_PATCH ${LIBSQUID_VERSION_STRING})
    message(STATUS "found version string ${LIBSQUID_VERSION_STRING}")

    if (NOT LIBSQUID_FIND_QUIETLY)
      message(STATUS "Found LIBSQUID ${LIBSQUID_VERSION_MAJOR}.${LIBSQUID_VERSION_MINOR}.${LIBSQUID_VERSION_PATCH}: ${LIBSQUID_LIBRARIES}")
    endif (NOT LIBSQUID_FIND_QUIETLY)
  else (LIBSQUID_FOUND)
    if (LIBSQUID_FIND_REQUIRED)
      message(STATUS "LIBSQUID not found.")
    endif (LIBSQUID_FIND_REQUIRED)
  endif (LIBSQUID_FOUND)

  mark_as_advanced(LIBSQUID_INCLUDE_DIR LIBSQUID_LIBRARIES)

endif (LIBSQUID_INCLUDE_DIR AND LIBSQUID_LIBRARIES)
