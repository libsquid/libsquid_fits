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

project(SQUIDFITS)
cmake_minimum_required(VERSION 2.8)

# Read squidfits.c and get version and release string
file(READ "${SQUIDFITS_SOURCE_DIR}/squidfits.c" SQUIDFITS_C)
string(FIND "${SQUIDFITS_C}" "SQUIDFITS_VERSION" SQUIDFITS_VERSION_CHECK)
if (SQUIDFITS_VERSION_CHECK EQUAL -1)
   message(FATAL_ERROR "version number not in squidfits.c")
endif()
string(REGEX REPLACE ".*#define SQUIDFITS_VERSION[^0-9]*([0-9]+)\\.([0-9]+)\\.([0-9]+).*" 
	     "\\1.\\2.\\3" SQUIDFITS_VERSION_STRING "${SQUIDFITS_C}")
if (SQUIDFITS_VERSION_STRING STREQUAL ifdef)
   message(FATAL_ERROR "version string not set in squidfits.c")
endif()
message("SQUIDFITS_VERSION_STRING ${SQUIDFITS_VERSION_STRING}")
string(FIND "${SQUIDFITS_C}" "SQUIDFITS_RELEASE" SQUIDFITS_RELEASE_CHECK)
if (SQUIDFITS_RELEASE_CHECK EQUAL -1)
   set(SQUIDFITS_RELEASE_STRING "")
else()
   string(REGEX REPLACE ".*#define SQUIDFITS_RELEASE[^0-9a-zA-Z\\.]*([0-9a-zA-Z\\.\\-]+).*" 
   	     "\\1" SQUIDFITS_RELEASE_STRING "${SQUIDFITS_C}")
   if (SQUIDFITS_RELEASE_STRING STREQUAL ifdef)
      set(SQUIDFITS_RELEASE_STRING "")
   endif()
   if (SQUIDFITS_RELEASE_STRING STREQUAL int)
      set(SQUIDFITS_RELEASE_STRING "")
   endif()
   message("SQUIDFITS_RELEASE_STRING ${SQUIDFITS_RELEASE_STRING}")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
		      "${CMAKE_SOURCE_DIR}/cmake_modules")

# Get list of source and header files
set(CMAKE_INCLUDE_CURRENT_DIR ON)
file(GLOB SQUIDFITS_SOURCES *.c)
file(GLOB SQUIDFITS_HEADERS *.h)

# Find necessary libraries
find_package(LIBSQUID REQUIRED)
message("LIBSQUID_FOUND = ${LIBSQUID_FOUND}")
message("LIBSQUID_INCLUDE_DIR = ${LIBSQUID_INCLUDE_DIR}")
message("LIBSQUID_LIBRARIES = ${LIBSQUID_LIBRARIES}")
if (NOT LIBSQUID_FOUND)
   message(FATAL_ERROR "LIBSQUID library not found")
endif() 
if (${LIBSQUID_VERSION_STRING} VERSION_LESS "0.5.3")
   message(FATAL_ERROR "LIBSQUID version ${LIBSQUID_VERSION_STRING} is too old")
endif()
find_package(LIBSQUIDWCS REQUIRED)
message("LIBSQUIDWCS_FOUND = ${LIBSQUIDWCS_FOUND}")
message("LIBSQUIDWCS_INCLUDE_DIR = ${LIBSQUIDWCS_INCLUDE_DIR}")
message("LIBSQUIDWCS_LIBRARIES = ${LIBSQUIDWCS_LIBRARIES}")
if (NOT LIBSQUIDWCS_FOUND)
   message(FATAL_ERROR "LIBSQUIDWCS library not found")
endif() 
if (${LIBSQUIDWCS_VERSION_STRING} VERSION_LESS "0.6")
   message(FATAL_ERROR "LIBSQUID_WCS version ${LIBSQUIDWCS_VERSION_STRING} is too old")
endif()
find_package(CFITSIO REQUIRED)
if (${CFITSIO_VERSION_STRING} VERSION_LESS "3.0")
   message(FATAL_ERROR "CFITSIO version ${CFITSIO_VERSION_STRING} is too old")
endif()
find_package(WCSLIB REQUIRED)
if (${WCSLIB_VERSION_STRING} VERSION_LESS "4.19")
   message(FATAL_ERROR "WCSLIB version ${WCSLIB_VERSION_STRING} is too old")
endif()
include_directories(${LIBSQUID_SOURCE_DIR}
	${LIBSQUID_INCLUDE_DIR}
	${LIBSQUIDWCS_SOURCE_DIR}
	${LIBSQUIDWCS_INCLUDE_DIR}
	${CFITSIO_INCLUDE_DIR}
	${WCSLIB_INCLUDE_DIR})
set(LIBS ${LIBS}
	 ${LIBSQUID_LIBRARIES}
	 ${LIBSQUIDWCS_LIBRARIES}
	 ${CFITSIO_LIBRARIES}
	 ${WCSLIB_LIBRARIES})
message("LIBS = ${LIBS}")
message("link dirs = ${LIBSQUIDWCS_BINARY_DIR} ${LIBSQUID_BINARY_DIR}")
link_directories(${LIBSQUIDWCS_BINARY_DIR} ${LIBSQUID_BINARY_DIR})
link_libraries(${LIBS} m)

foreach(infile ${SQUIDFITS_SOURCES})
   string(REGEX MATCH "[^/\\]+\\.c$" outfile "${infile}")
   string(REPLACE ".c" "" outfile ${outfile})
   message("adding target ${outfile} ${infile}")
   set(SQUIDFITS_TARGETS "${SQUIDFITS_TARGETS} ${outfile}")
   add_executable(${outfile} ${infile})
   install(TARGETS ${outfile}
	   DESTINATION bin
	   COMPONENT squidfits
   )
endforeach(infile)

set(CPACK_PACKAGE_NAME "squidfits")
set(CPACK_PACKAGE_VERSION ${SQUIDFITS_VERSION_STRING})
if (SQUIDFITS_RELEASE_STRING STREQUAL "")
   set(CPACK_PACKAGE_RELEASE "")
else()
   set(CPACK_PACKAGE_RELEASE "~${SQUIDFITS_RELEASE_STRING}")
endif()
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY 
    "Programs to use LibSQUID library to manipulate FITS images")
set(CPACK_PACKAGE_CONTACT "Jim Wren; jwren@lanl.gov")
set(CPACK_PACKAGE_VENDOR ${CPACK_PACKAGE_CONTACT})
if (NOT SQUIDFITS_RELEASE_STRING STREQUAL "")
   set(CPACK_RPM_PACKAGE_RELEASE "${SQUIDFITS_RELEASE_STRING}")
endif()
set(CPACK_RPM_PACKAGE_REQUIRES
   "libsquid >= 0.5.3, libsquid_wcs > 0.6, cfitsio >= 3.0, wcslib >= 4.19")
set(CPACK_DEBIAN_PACKAGE_DEPENDS
   "libsquid (>= 0.5.3), libsquid_wcs (> 0.6) cfitsio (>= 3.0), wcslib (>= 4.19)")
set(CPACK_GENERATOR "RPM" "DEB")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "~$;[.]bak$/;[.]swp$;[.]svn$")
set(CPACK_PACKAGE_DIRECTORY ${CMAKE_BINARY_DIR}/packages)
set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}${CPACK_PACKAGE_RELEASE}.${CMAKE_SYSTEM_PROCESSOR}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}${CPACK_PACKAGE_RELEASE}")
include(CPack)
