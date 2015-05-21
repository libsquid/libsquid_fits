# - Try to find CFITSIO
# Once done this will define
#
#  CFITSIO_FOUND - system has CFITSIO
#  CFITSIO_INCLUDEDIR - the CFITSIO include directory
#  CFITSIO_LIBRARIES - Link these to use CFITSIO
#  CFITSIO_VERSION_STRING - Human readable version number of cfitsio
#  CFITSIO_VERSION_MAJOR  - Major version number of cfitsio
#  CFITSIO_VERSION_MINOR  - Minor version number of cfitsio

# Copyright (c) 2006, Jasem Mutlaq <mutlaqja@ikarustech.com>
# Based on FindLibfacile by Carsten Niehaus, <cniehaus@gmx.de>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

# Force CFITSIO_LIBDIR availability in ccmake
set(CFITSIO_LIBDIR ${CFITSIO_LIBDIR} CACHE STRING "CFITSIO library path")

mark_as_advanced(CFITSIO_INCLUDEDIR CFITSIO_LIBDIR CFITSIO_LIBRARIES)

# Assume not found
set(CFITSIO_FOUND FALSE)

# Call pkg-config
find_package(PkgConfig)
find_package(PkgConfig)
if (PKG_CONFIG_FOUND)
  if (UNIX)
    pkg_check_modules(PC_CFITSIO QUIET cfitsio)
  endif()
endif()

# Form CFITSIO_INCLUDEDIR
if (CFITSIO_INCLUDEDIR)
  find_path(CFITSIO_INCLUDEDIR fitsio.h
            # JM: Packages from different distributions have different suffixes
            PATH_SUFFIXES libcfitsio3 libcfitsio0 cfitsio
            PATHS ${CFITSIO_INCLUDEDIR})
else()
  find_path(CFITSIO_INCLUDEDIR fitsio.h
            # JM: Packages from different distributions have different suffixes
            PATH_SUFFIXES libcfitsio3 libcfitsio0 cfitsio
            HINTS ${PC_CFITSIO_INCLUDEDIR} ${PC_CFITSIO_INCLUDE_DIRS}
            PATHS ${CFITSIO_INCLUDEDIR} $ENV{CFITSIO_INCLUDEDIR})
endif()
if (CFITSIO_INCLUDEDIR STREQUAL "CFITSIO_INCLUDEDIR-NOTFOUND")
  if (CFITSIO_FIND_REQUIRED)
    message(FATAL_ERROR "CFITSIO_INCLUDEDIR not found.")
  else()
    if (NOT CFITSIO_FIND_QUIETLY)
      message(STATUS "CFITSIO_INCLUDEDIR not found.")
    endif()
    return()
  endif()
endif()

# Form CFITSIO_VERSION_STRING
FILE(STRINGS "${CFITSIO_INCLUDEDIR}/fitsio.h" CFITSIO_VERSION_STRING
     REGEX "^[ \t]*#define[ \t]+CFITSIO_VERSION[ \t]+[0-9]+\\.[0-9]+[ \t]*$")
if (CFITSIO_VERSION_STRING STREQUAL "")
  if (CFITSIO_FIND_REQUIRED)
    message(FATAL_ERROR "CFITSIO_VERSION_STRING not found.")
  else()
    if (NOT CFITSIO_FIND_QUIETLY)
      message(STATUS "CFITSIO_VERSION_STRING not found.")
    endif()
    return()
  endif()
endif()
STRING(REGEX REPLACE ".*CFITSIO_VERSION[ \t]+([0-9]+)\\.([0-9]+).*" "\\1.\\2" CFITSIO_VERSION_STRING "${CFITSIO_VERSION_STRING}")
STRING(REGEX REPLACE "([0-9]+)\\.([0-9]+)" "\\1" CFITSIO_VERSION_MAJOR ${CFITSIO_VERSION_STRING})
STRING(REGEX REPLACE "([0-9]+)\\.([0-9]+)" "\\2" CFITSIO_VERSION_MINOR ${CFITSIO_VERSION_STRING})

# Form CFITSIO_LIBRARIES based on CFITSIO_LIBDIR
if (CFITSIO_LIBDIR)
  find_library(CFITSIO_LIBRARIES NAMES cfitsio
               PATHS ${CFITSIO_LIBDIR})
else()
  find_library(CFITSIO_LIBRARIES NAMES cfitsio
               HINTS ${PC_CFITSIO_LIBDIR} ${PC_CFITSIO_LIBRARY_DIRS}
               PATHS ${CFITSIO_LIBDIR} $ENV{CFITSIO_LIBDIR})
endif()
if (CFITSIO_LIBRARIES STREQUAL "CFITSIO_LIBRARIES-NOTFOUND")
  if ((CFITSIO_LIBDIR STREQUAL "") AND (NOT CFITSIO_FIND_QUIETLY))
      message(WARNING "CFITSIO_LIBDIR empty")
  endif()
  if (CFITSIO_FIND_REQUIRED)
    message(FATAL_ERROR "CFITSIO_LIBRARIES not found.")
  else()
    if (NOT CFITSIO_FIND_QUIETLY)
      message(STATUS "CFITSIO_LIBRARIES not found.")
    endif()
    return()
  endif()
endif()

if (NOT CFITSIO_FIND_QUIETLY)
  message(STATUS "Found CFITSIO ${CFITSIO_VERSION_MAJOR}.${CFITSIO_VERSION_MINOR}: ${CFITSIO_LIBRARIES}")
endif()

set(CFITSIO_FOUND TRUE)
