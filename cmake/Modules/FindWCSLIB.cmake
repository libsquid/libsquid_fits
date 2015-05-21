# - Try to find WCSLIB
# Once done this will define
#
#  WCSLIB_FOUND - system has WCSLIB
#  WCSLIB_INCLUDEDIR - the WCSLIB include directory
#  WCSLIB_LIBRARIES - Link these to use WCSLIB
#  WCSLIB_VERSION_STRING - Human readable version number
#  WCSLIB_VERSION_MAJOR - major version number
#  WCSLIB_VERSION_MINOR - minor version number

# Copyright (c) 2006, Jasem Mutlaq <mutlaqja@ikarustech.com>
# Based on FindLibfacile by Carsten Niehaus, <cniehaus@gmx.de>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

# Force WCSLIB_LIBDIR availability in ccmake
set(WCSLIB_LIBDIR ${WCSLIB_LIBDIR} CACHE STRING "WCSLIB library path")

mark_as_advanced(WCSLIB_INCLUDEDIR WCSLIB_LIBDIR WCSLIB_LIBRARIES)

# Assume not found
set(WCSLIB_FOUND FALSE)

# Call pkg-config
find_package(PkgConfig)
if (PKG_CONFIG_FOUND)
  if (UNIX)
    pkg_check_modules(PC_WCSLIB QUIET wcslib)
  endif()
endif()

# Form WCSLIB_INCLUDEDIR
if (WCSLIB_INCLUDEDIR)
  find_path(WCSLIB_INCLUDEDIR wcs.h
            PATH_SUFFIXES wcslib
            PATHS ${WCSLIB_INCLUDEDIR})
else()
  find_path(WCSLIB_INCLUDEDIR wcs.h
            PATH_SUFFIXES wcslib
            HINTS ${PC_WCSLIB_INCLUDEDIR} ${PC_WCSLIB_INCLUDE_DIRS}
            PATHS ${WCSLIB_INCLUDEDIR} $ENV{WCSLIB_INCLUDEDIR})
endif()
if (WCSLIB_INCLUDEDIR STREQUAL "WCSLIB_INCLUDEDIR-NOTFOUND")
  if (WCSLIB_FIND_REQUIRED)
    message(FATAL_ERROR "WCSLIB_INCLUDEDIR not found.")
  else()
    if (NOT WCSLIB_FIND_QUIETLY)
      message(STATUS "WCSLIB_INCLUDEDIR not found.")
    endif()
    return()
  endif()
endif()

# Form WCSLIB_VERSION_STRING
FILE(STRINGS "${WCSLIB_INCLUDEDIR}/wcs.h" WCSLIB_VERSION_STRING
     REGEX "^.*wcs.h,v[ \t]+[0-9]+\\.[0-9]+[ \t]+.*$")
if (WCSLIB_VERSION_STRING STREQUAL "")
  if (WCSLIB_FIND_REQUIRED)
    message(FATAL_ERROR "WCSLIB_VERSION_STRING not found.")
  else()
    if (NOT WCSLIB_FIND_QUIETLY)
      message(STATUS "WCSLIB_VERSION_STRING not found.")
    endif()
    return()
  endif()
endif()
STRING(REGEX REPLACE ".*wcs.h,v[ \t]+([0-9]+)\\.([0-9]+).*" "\\1.\\2" WCSLIB_VERSION_STRING "${WCSLIB_VERSION_STRING}")
STRING(REGEX REPLACE "([0-9]+)\\.([0-9]+)" "\\1" WCSLIB_VERSION_MAJOR ${WCSLIB_VERSION_STRING})
STRING(REGEX REPLACE "([0-9]+)\\.([0-9]+)" "\\2" WCSLIB_VERSION_MINOR ${WCSLIB_VERSION_STRING})

# Form WCSLIB_LIBRARIES based on WCSLIB_LIBDIR
if (WCSLIB_LIBDIR)
  find_library(WCSLIB_LIBRARIES NAMES wcs
               PATHS ${WCSLIB_LIBDIR})
else()
  find_library(WCSLIB_LIBRARIES NAMES wcs
               HINTS ${PC_WCSLIB_LIBDIR} ${PC_WCSLIB_LIBRARY_DIRS}
               PATHS ${WCSLIB_LIBDIR} $ENV{WCSLIB_LIBDIR})
endif()
if (WCSLIB_LIBRARIES STREQUAL "WCSLIB_LIBRARIES-NOTFOUND")
  if ((WCSLIB_LIBDIR STREQUAL "") AND (NOT WCSLIB_FIND_QUIETLY))
      message(WARNING "WCSLIB_LIBDIR empty")
  endif()
  if (WCSLIB_FIND_REQUIRED)
    message(FATAL_ERROR "WCSLIB_LIBRARIES not found.")
  else()
    if (NOT WCSLIB_FIND_QUIETLY)
      message(STATUS "WCSLIB_LIBRARIES not found.")
    endif()
    return()
  endif()
endif()

if (NOT WCSLIB_FIND_QUIETLY)
  message(STATUS "Found WCSLIB ${WCSLIB_VERSION_MAJOR}.${WCSLIB_VERSION_MINOR}: ${WCSLIB_LIBRARIES}")
endif()

set(WCSLIB_FOUND TRUE)
