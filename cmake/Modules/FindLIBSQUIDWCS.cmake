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

# Force LIBSQUIDWCS_LIBDIR availability in ccmake
set(LIBSQUIDWCS_LIBDIR ${LIBSQUIDWCS_LIBDIR} CACHE STRING "LibSQUID_WCS library path")

mark_as_advanced(LIBSQUIDWCS_INCLUDEDIR LIBSQUIDWCS_LIBDIR LIBSQUIDWCS_LIBRARIES)

# Assume not found
set(LIBSQUIDWCS_FOUND FALSE)

# Call pkg-config
find_package(PkgConfig)
if (PKG_CONFIG_FOUND)
  if (UNIX)
    pkg_check_modules(PC_LIBSQUID_WCS QUIET libsquid_wcs)
  endif()
endif()

# Form LIBSQUIDWCS_INCLUDEDIR
if (LIBSQUIDWCS_INCLUDEDIR)
  find_path(LIBSQUIDWCS_INCLUDEDIR libsquid_wcs.h
            PATHS ${LIBSQUIDWCS_INCLUDEDIR})
else()
  find_path(LIBSQUIDWCS_INCLUDEDIR libsquid_wcs.h
            HINTS ${PC_LIBSQUIDWCS_INCLUDEDIR} ${PC_LIBSQUIDWCS_INCLUDE_DIRS}
            PATHS ${LIBSQUIDWCS_INCLUDEDIR} $ENV{LIBSQUIDWCS_INCLUDEDIR})
endif()
if (LIBSQUIDWCS_INCLUDEDIR STREQUAL "LIBSQUIDWCS_INCLUDEDIR-NOTFOUND")
  if (LIBSQUIDWCS_FIND_REQUIRED)
    message(FATAL_ERROR "LIBSQUIDWCS_INCLUDEDIR not found.")
  else()
    if (NOT LIBSQUIDWCS_FIND_QUIETLY)
      message(STATUS "LIBSQUIDWCS_INCLUDEDIR not found.")
    endif()
    return()
  endif()
endif()

# Form LIBSQUIDWCS_VERSION_STRING
FILE(STRINGS "${LIBSQUIDWCS_INCLUDEDIR}/libsquid_wcs.h" LIBSQUIDWCS_VERSION_STRING
     REGEX "^[ \t]*#define[ \t]+LIBSQUIDWCS_VERSION[ \t]+\"[0-9]+\\.[0-9]+\\.[0-9]+\"[ \t]*$")
if (LIBSQUIDWCS_VERSION_STRING STREQUAL "")
  if (LIBSQUIDWCS_FIND_REQUIRED)
    message(FATAL_ERROR "LIBSQUIDWCS_VERSION_STRING not found.")
  else()
    if (NOT LIBSQUIDWCS_FIND_QUIETLY)
      message(STATUS "LIBSQUIDWCS_VERSION_STRING not found.")
    endif()
    return()
  endif()
endif()
STRING(REGEX REPLACE ".*\"([0-9]+)\\.([0-9]+)\\.([0-9]+)\".*" "\\1.\\2.\\3" LIBSQUIDWCS_VERSION_STRING "${LIBSQUIDWCS_VERSION_STRING}")
STRING(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\1" LIBSQUIDWCS_VERSION_MAJOR ${LIBSQUIDWCS_VERSION_STRING})
STRING(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\2" LIBSQUIDWCS_VERSION_MINOR ${LIBSQUIDWCS_VERSION_STRING})
STRING(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\3" LIBSQUIDWCS_VERSION_PATCH ${LIBSQUIDWCS_VERSION_STRING})

# Form LIBSQUIDWCS_LIBRARIES based on LIBSQUIDWCS_LIBDIR
if (LIBSQUIDWCS_LIBDIR)
  find_library(LIBSQUIDWCS_LIBRARIES NAMES squid_wcs
               PATHS ${LIBSQUIDWCS_LIBDIR})
else()
  find_library(LIBSQUIDWCS_LIBRARIES NAMES squid_wcs
               HINTS ${PC_LIBSQUIDWCS_LIBDIR} ${PC_LIBSQUIDWCS_LIBRARY_DIRS}
               PATHS ${LIBSQUIDWCS_LIBDIR} $ENV{LIBSQUIDWCS_LIBDIR})
endif()
if (LIBSQUIDWCS_LIBRARIES STREQUAL "LIBSQUIDWCS_LIBRARIES-NOTFOUND")
  if ((LIBSQUIDWCS_LIBDIR STREQUAL "") AND (NOT LIBSQUIDWCS_FIND_QUIETLY))
      message(WARNING "LIBSQUIDWCS_LIBDIR empty")
  endif()
  if (LIBSQUIDWCS_FIND_REQUIRED)
    message(FATAL_ERROR "LIBSQUIDWCS_LIBRARIES not found.")
  else()
    if (NOT LIBSQUIDWCS_FIND_QUIETLY)
      message(STATUS "LIBSQUIDWCS_LIBRARIES not found.")
    endif()
    return()
  endif()
endif()

if (NOT LIBSQUIDWCS_FIND_QUIETLY)
  message(STATUS "Found LIBSQUIDWCS ${LIBSQUIDWCS_VERSION_MAJOR}.${LIBSQUIDWCS_VERSION_MINOR}.${LIBSQUIDWCS_VERSION_PATCH}: ${LIBSQUIDWCS_LIBRARIES}")
endif()

set(LIBSQUIDWCS_FOUND TRUE)
