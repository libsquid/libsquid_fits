# - Try to find LIBSQUID
# Once done this will define
#
#  LIBSQUID_FOUND - system has LIBSQUID
#  LIBSQUID_INCLUDEDIR - the LIBSQUID include directory
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

# Force LIBSQUID_LIBDIR availability in ccmake
set(LIBSQUID_LIBDIR ${LIBSQUID_LIBDIR} CACHE STRING "LibSQUID library path")

mark_as_advanced(LIBSQUID_INCLUDEDIR LIBSQUID_LIBDIR LIBSQUID_LIBRARIES)

# Assume not found
set(LIBSQUID_FOUND FALSE)

# Call pkg-config
find_package(PkgConfig)
if (PKG_CONFIG_FOUND)
  if (UNIX)
    pkg_check_modules(PC_LIBSQUID QUIET libsquid)
  endif()
endif()

# Form LIBSQUID_INCLUDEDIR
if (LIBSQUID_INCLUDEDIR)
  find_path(LIBSQUID_INCLUDEDIR libsquid.h
            PATHS ${LIBSQUID_INCLUDEDIR})
else()
  find_path(LIBSQUID_INCLUDEDIR libsquid.h
            HINTS ${PC_LIBSQUID_INCLUDEDIR} ${PC_LIBSQUID_INCLUDE_DIRS}
            PATHS ${LIBSQUID_INCLUDEDIR} $ENV{LIBSQUID_INCLUDEDIR})
endif()
if (LIBSQUID_INCLUDEDIR STREQUAL "LIBSQUID_INCLUDEDIR-NOTFOUND")
  if (LIBSQUID_FIND_REQUIRED)
    message(FATAL_ERROR "LIBSQUID_INCLUDEDIR not found.")
  else()
    if (NOT LIBSQUID_FIND_QUIETLY)
      message(STATUS "LIBSQUID_INCLUDEDIR not found.")
    endif()
    return()
  endif()
endif()

# Form LIBSQUID_VERSION_STRING
FILE(STRINGS "${LIBSQUID_INCLUDEDIR}/libsquid.h" LIBSQUID_VERSION_STRING
     REGEX "^[ \t]*#define[ \t]+LIBSQUID_VERSION[ \t]+\"[0-9]+\\.[0-9]+\\.[0-9]+\"[ \t]*$")
if (LIBSQUID_VERSION_STRING STREQUAL "")
  if (LIBSQUID_FIND_REQUIRED)
    message(FATAL_ERROR "LIBSQUID_VERSION_STRING not found.")
  else()
    if (NOT LIBSQUID_FIND_QUIETLY)
      message(STATUS "LIBSQUID_VERSION_STRING not found.")
    endif()
    return()
  endif()
endif()
STRING(REGEX REPLACE ".*\"([0-9]+)\\.([0-9]+)\\.([0-9]+)\".*" "\\1.\\2.\\3" LIBSQUID_VERSION_STRING "${LIBSQUID_VERSION_STRING}")
STRING(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\1" LIBSQUID_VERSION_MAJOR ${LIBSQUID_VERSION_STRING})
STRING(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\2" LIBSQUID_VERSION_MINOR ${LIBSQUID_VERSION_STRING})
STRING(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\3" LIBSQUID_VERSION_PATCH ${LIBSQUID_VERSION_STRING})

# Form LIBSQUID_LIBRARIES based on LIBSQUID_LIBDIR
if (LIBSQUID_LIBDIR)
  find_library(LIBSQUID_LIBRARIES NAMES squid
               PATHS ${LIBSQUID_LIBDIR})
else()
  find_library(LIBSQUID_LIBRARIES NAMES squid
               HINTS ${PC_LIBSQUID_LIBDIR} ${PC_LIBSQUID_LIBRARY_DIRS}
               PATHS ${LIBSQUID_LIBDIR} $ENV{LIBSQUID_LIBDIR})
endif()
if (LIBSQUID_LIBRARIES STREQUAL "LIBSQUID_LIBRARIES-NOTFOUND")
  if ((LIBSQUID_LIBDIR STREQUAL "") AND (NOT LIBSQUID_FIND_QUIETLY))
      message(WARNING "LIBSQUID_LIBDIR empty")
  endif()
  if (LIBSQUID_FIND_REQUIRED)
    message(FATAL_ERROR "LIBSQUID_LIBRARIES not found.")
  else()
    if (NOT LIBSQUID_FIND_QUIETLY)
      message(STATUS "LIBSQUID_LIBRARIES not found.")
    endif()
    return()
  endif()
endif()

if (NOT LIBSQUID_FIND_QUIETLY)
  message(STATUS "Found LIBSQUID ${LIBSQUID_VERSION_MAJOR}.${LIBSQUID_VERSION_MINOR}.${LIBSQUID_VERSION_PATCH}: ${LIBSQUID_LIBRARIES}")
endif()

set(LIBSQUID_FOUND TRUE)
