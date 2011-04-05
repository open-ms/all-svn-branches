# -*- mode: C++; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

include(CppcheckTargets)

if(CPPCHECK_FOUND)
  set(SOURCE_FILE_REGEX "\\.C$")

  # library checks
  set(OpenMS_cppcheck_sources)
  foreach(i ${OpenMS_sources})
    string( REGEX MATCH ${SOURCE_FILE_REGEX} is_source_file ${i} )
    if(is_source_file)
      #add_cppcheck_sources(${i} ${i} STYLE FAIL_ON_WARNINGS)
      list(APPEND OpenMS_cppcheck_sources ${i})
    endif()
  endforeach()

  add_cppcheck_sources("OpenMS" ${OpenMS_cppcheck_sources} STYLE FAIL_ON_WARNINGS)

  # GUI library checks
  set(OpenMSVisual_cppcheck_sources)
  foreach(i ${OpenMSVisual_sources})
    string( REGEX MATCH ${SOURCE_FILE_REGEX} is_source_file ${i} )
    if(is_source_file)
      list(APPEND OpenMSVisual_cppcheck_sources ${i})
      # add_cppcheck_sources(${i} ${i} STYLE FAIL_ON_WARNINGS)
    endif()
  endforeach()

  add_cppcheck_sources("OpenMS_GUI" ${OpenMSVisual_cppcheck_sources} STYLE FAIL_ON_WARNINGS)

  # TOPP checks
  foreach(i ${TOPP_executables})
    add_cppcheck( ${i} STYLE FAIL_ON_WARNINGS )
  endforeach()

  # UTILS checks
  foreach(i ${UTILS_executables})
    add_cppcheck( ${i} STYLE FAIL_ON_WARNINGS )
  endforeach()

  # TEST checks


else()
  message(STATUS "Missing CPPCHECK executable .. Abort CppCheck")
endif()