# CMake installation file
# Created 2008/10 during OpenMS retreat by CB && AB

PROJECT("OpenMS")

########################################################
### entries meant to be configured using CMake cache ###
### - do NOT hardcode them here!										 ###
### - edit them within CMakeCache.txt								 ### 
########################################################
# CUSTOM_CONTRIB_DIR
# OPENMS_DEBUG
# OPENMS_ENABLE_EXTERNALMEMORY
# CF_DB_PLUGIN

########################################################
###    manual entries (edit this for new release)    ###
########################################################

set(CF_PACKAGE_VERSION "1.3" CACHE INTERNAL "OpenMS VERSION" FORCE)

########################################################
###    end manual																	   ###
########################################################

# Heart of the BUILD system : only edit when you know what you´re doing (we don´t)
# quick manual for most commands: http://www.cmake.org/cmake/help/cmake2.6docs.html

cmake_minimum_required(VERSION 2.6)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

## define some directories
if (CMAKE_INSTALL_PREFIX)
	set(DIR_SHARE ${CMAKE_INSTALL_PREFIX}/share/ CACHE INTERNAL "Path to the /share directory after installation.")
else (CMAKE_INSTALL_PREFIX)
	set(DIR_SHARE ${PROJECT_SOURCE_DIR}/share/ CACHE INTERNAL "Path to the /share directory after installation.")
endif (CMAKE_INSTALL_PREFIX)

########################################################
###    compiler flags																 ###
########################################################

## general
#TODO confirm that these flags are general
add_definitions("-Wall -pedantic -Wno-long-long -Wno-non-virtual-dtor")

## platform dependent compiler flags:
include(CheckCXXCompilerFlag)
if (NOT WIN32) # we only want fPIC on non-windows systems (fPIC is implicitly true there)
	CHECK_CXX_COMPILER_FLAG("-fPIC" WITH_FPIC) 
	if (WITH_FPIC)
		add_definitions(-fPIC)
	endif()
endif()


########################################################
###    external libs																 ###
########################################################

## contrib

## disable AndiMS for 64-bit platforms (until fixed)
if (CMAKE_SIZEOF_VOID_P MATCHES "8")
	set(USE_ANDIMS "false" CACHE INTERNAL "Support AndiMS?!")
else()
	set(USE_ANDIMS "true" CACHE INTERNAL "Support AndiMS?!")
endif()

set(CUSTOM_CONTRIB_DIR CACHE PATH "User defined location of contrib dir. If left empty we assume the contrib to be in OpenMS/contrib!")
set(CONTRIB_DIR ${PROJECT_SOURCE_DIR}/contrib/ CACHE INTERNAL "Final contrib path after looking at custom_contrib_path. defaults to OpenMS/contrib")

#message("CUSTOM contrib is: " ${CUSTOM_CONTRIB_DIR}) 
 
if (CUSTOM_CONTRIB_DIR STREQUAL "")
	message(STATUS "Using DEFAULT setting for contrib directory: ${CONTRIB_DIR}")
else()
	set(CONTRIB_DIR ${CUSTOM_CONTRIB_DIR})
	message(STATUS "Using user defined contrib directory: ${CONTRIB_DIR}")
endif()

FIND_PATH(CONFIRM_CONTRIB_DIR README.txt PATHS ${CONTRIB_DIR} NO_DEFAULT_PATH)
if (CONFIRM_CONTRIB_DIR STREQUAL "CONFIRM_CONTRIB_DIR-NOTFOUND")
	message(FATAL_ERROR "The path to the Contrib directory (${CONTRIB_DIR}) is not correct!")
endif()

set(CONTRIB_INCLUDE_DIR ${CONTRIB_DIR}/include/ CACHE INTERNAL "contrib include dir")
set(CONTRIB_LIB_DIR ${CONTRIB_DIR}/lib/ CACHE INTERNAL "contrib lib dir")


# TODO check whether this is necessary
FIND_PATH(GSL_INCLUDE_DIR gsl/gsl_cdf.h NAMES gsl/gsl_randist.h PATHS ${CONTRIB_INCLUDE_DIR} DOC "GSL header include dir")
if (GSL_INCLUDE_DIR STREQUAL "GSL_INCLUDE_DIR-NOTFOUND")
	MESSAGE(FATAL_ERROR "Unable to find GSL include dir! Please make sure it is part of the contrib (which we assume to be at: ${CONTRIB_DIR}")
else()
	MESSAGE(STATUS "Found GSL include dir at: " ${GSL_INCLUDE_DIR})
endif()

# TODO check whether this is necessary
#FIND_LIBRARY(GSL_LIBRARY NAMES acctres hal hal.dll PATHS "c:/windows/system32/" )
#FIND_LIBRARY(GSL_LIBRARY NAMES gsl libgsl libgsl.dll PATHS  ${CONTRIB_LIB_DIR} DOC "GSL library dir")
if (GSL_LIBRARY STREQUAL "GSL_LIBRARY-NOTFOUND")
	MESSAGE(FATAL_ERROR "Unable to find GSL library! Please make sure it is part of the contrib (which we assume to be at: ${CONTRIB_DIR}")
else()
	MESSAGE(STATUS "Found GSL library at: " ${GSL_LIBRARY})
endif()

# TODO check for gslcblaslib

#FIND_LIBRARY(CGAL_LIBRARY NAMES CGAL libCGAL PATHS  ${CONTRIB_LIB_DIR} DOC "CGAL library dir")
if (CGAL_LIBRARY STREQUAL "CGAL_LIBRARY-NOTFOUND")
	MESSAGE(FATAL_ERROR "Unable to find CGAL library! Please make sure it is part of the contrib (which we assume to be at: ${CONTRIB_DIR}")
else()
	MESSAGE(STATUS "Found CGAL library at: " ${CGAL_LIBRARY})
endif()

#FIND_LIBRARY(XERCESC_LIBRARY NAMES xerces-c libxerces-c PATHS ${CONTRIB_LIB_DIR} DOC "xerces-c library dir")
if (XERCESC_LIBRARY STREQUAL "XERCESC_LIBRARY-NOTFOUND")
	MESSAGE(FATAL_ERROR "Unable to find xerces-c library! Please make sure it is part of the contrib (which we assume to be at: ${CONTRIB_DIR}")
else()
	MESSAGE(STATUS "Found xerces-c library at: " ${XERCESC_LIBRARY})
endif()

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include/)
INCLUDE_DIRECTORIES(${CONTRIB_INCLUDE_DIR})

## QT
SET(QT_USE_QTXML 1)
SET(QT_USE_QTNETWORK 1)
SET(QT_USE_QTSQL 1)
SET(QT_USE_QTOPENGL 1)
FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})
INCLUDE(UseQt4)

########################################################
###    useful programms															 ###
########################################################

find_program(PROG_Svnversion             
             svnversion
             PATHS "c:/programme/subversion/bin" "c:/program files/subversion/bin"   ### additional search paths (along with $PATH)
             DOC "svnversion executable which helps in determining the svn revision when building TOPP tools")

if (PROG_Svnversion STREQUAL "PROG_Svnversion-NOTFOUND")
	message(STATUS "The programm svnversion could not be found. SVN-revision information will not be available! Add the location of svnversion(.exe) to your PATH environment variable if you require SVN-revision.")
endif()

########################################################
###    configure config.h														 ###
########################################################


# user DEBUG flag (needs a cmake flag/cache entry, so we define it empty which causes it to show up in the wizard)
set(OPENMS_DEBUG CACHE BOOL "Enable debug output for certain classes")

# can be set automatically
set(CF_OPENMS_DATA_PATH ${DIR_SHARE}OpenMS-${CF_PACKAGE_VERSION}/ CACHE INTERNAL "Path to the shared documents of OpenMS.")

# check for Microsoft Visual Studio compiler
if (MSVC)
	set(OPENMS_COMPILER_MSVC "1" CACHE INTERNAL "Do we use Microsoft Compiler?")
endif()
# check for G++
if (CMAKE_COMPILER_IS_GNUCXX)
	set(OPENMS_COMPILER_GXX "1" CACHE INTERNAL "Do we use G++ Compiler?")
endif()

INCLUDE(TestBigEndian)
TEST_BIG_ENDIAN(OPENMS_BIG_ENDIAN)

# check 32/64 bit architecture
if (CMAKE_SIZEOF_VOID_P MATCHES "8")
	set(OPENMS_64BIT_ARCHITECTURE, "1")
endif()


include(CheckTypeSize) #: Check sizeof a type
CHECK_TYPE_SIZE("unsigned char" SIZE_UCHAR)
CHECK_TYPE_SIZE("unsigned short" SIZE_USHORT)
CHECK_TYPE_SIZE("unsigned int" SIZE_UINT)
CHECK_TYPE_SIZE("unsigned long" SIZE_ULONG)
CHECK_TYPE_SIZE("unsigned long long" SIZE_ULONGLONG)
CHECK_TYPE_SIZE("short" SIZE_SHORT)
CHECK_TYPE_SIZE("int" SIZE_INT)
CHECK_TYPE_SIZE("long" SIZE_LONG)
CHECK_TYPE_SIZE("long long" SIZE_LONGLONG)


CHECK_TYPE_SIZE("int32_t" SIZE_INT32)
if (HAVE_SIZE_INT32)
	set(CF_OPENMS_INT32_TYPE int32_t)
else()
	# search for another Int32 type
	if (SIZE_INT MATCHES "4")
		set(CF_OPENMS_INT32_TYPE int)
	elseif (SIZE_SHORT MATCHES "4")
		set(CF_OPENMS_INT32_TYPE short)	
	elseif (SIZE_LONG MATCHES "4")
		set(CF_OPENMS_INT32_TYPE long)	
	else()	
		message(FATAL_ERROR "Cannot find signed 32bit integer type. Please contact the developers!")		
	endif()
endif()

CHECK_TYPE_SIZE("int64_t" SIZE_INT64)
if (HAVE_SIZE_INT64)
	set(CF_OPENMS_INT64_TYPE int64_t)
else()
	# search for another Int64 type
	if (SIZE_INT MATCHES "8")
		set(CF_OPENMS_INT64_TYPE int)
	elseif (SIZE_LONG MATCHES "8")
		set(CF_OPENMS_INT64_TYPE long)	
	elseif (SIZE_LONGLONG MATCHES "8")
		set(CF_OPENMS_INT64_TYPE "long long")	
	else()
		message(FATAL_ERROR "Cannot find signed 64bit integer type. Please contact the developers!")		
	endif()
endif()

CHECK_TYPE_SIZE("uint8_t" SIZE_UINT8)
if (HAVE_SIZE_UINT8)
	set(CF_OPENMS_BYTE_TYPE uint8_t)
else()
	# search for another uint8 type
	if (SIZE_UCHAR MATCHES "1")
		set(CF_OPENMS_BYTE_TYPE "unsigned char")
	elseif (SIZE_USHORT MATCHES "1")
		set(CF_OPENMS_BYTE_TYPE "unsigned short")	
	else()
		message(FATAL_ERROR "Cannot find unsigned 8bit integer (byte) type. Please contact the developers!")		
	endif()
endif()


CHECK_TYPE_SIZE("uint64_t" SIZE_UINT64)
if (HAVE_SIZE_UINT64)
	set(CF_OPENMS_UINT64_TYPE uint64_t)
else()
	# search for another uint64 type
	if (SIZE_ULONG MATCHES "8")
		set(CF_OPENMS_UINT64_TYPE "unsigned long")
	elseif (SIZE_ULONGLONG MATCHES "8")
		set(CF_OPENMS_UINT64_TYPE "unsigned long long")	
	else()
		message(FATAL_ERROR "Cannot find uint64 type. Please contact the developers!")		
	endif()
endif()

## system headers:
include(CheckIncludeFileCXX) # Check if the include file exists.

CHECK_INCLUDE_FILE_CXX("unistd.h" OPENMS_HAS_UNISTD_H)
CHECK_INCLUDE_FILE_CXX("process.h" OPENMS_HAS_PROCESS_H)

CHECK_INCLUDE_FILE_CXX("time.h" OPENMS_HAS_TIME_H)
CHECK_INCLUDE_FILE_CXX("sys/types.h" OPENMS_HAS_SYS_TYPES_H)
CHECK_INCLUDE_FILE_CXX("sys/times.h" OPENMS_HAS_SYS_TIMES_H)
CHECK_INCLUDE_FILE_CXX("sys/time.h"  OPENMS_HAS_SYS_TIME_H)

# user flag (needs a cmake flag/cache entry, so we define it empty which causes it to show up in the wizard)
set(OPENMS_ENABLE_EXTERNALMEMORY CACHE BOOL "User switch to enable external memory support (if you run out of RAM+swap). Warning: performance penalty!")

include(CheckFunctionExists)
## in MinGW we have the signal.h header, but no kill() as in Linux, so we need to check for the kill() function 
CHECK_FUNCTION_EXISTS("kill" OPENMS_HAS_KILL)
CHECK_FUNCTION_EXISTS("sysconf" OPENMS_HAS_SYSCONF)

## ANDIMS flag USE_ANDIMS defined above!
if (NOT DEFINED USE_ANDIMS)
	message(FATAL_ERROR "Cmake script was re-ordered and is now invalid! Please make sure that USE_ANIMS is defined when config.h.in is configured!")		
endif()

# user flag with default "QMYSQL" (needs a cmake flag/cache entry, so we define it empty which causes it to show up in the wizard)
set(CF_DB_PLUGIN "QMYSQL" CACHE BOOL "User switch to change the DB plugin.")


## replace any variables in config.h.in with current values
configure_file(${PROJECT_SOURCE_DIR}/source/config/config.h.in ${PROJECT_SOURCE_DIR}/include/OpenMS/config.h)

########################################################
###    configure Doxyfile														 ###
########################################################

## CF_PACKAGE_VERSION (already defined above in manual section)

set(CF_OPENMS_PATH ${PROJECT_BINARY_DIR})

## replace any variables in config.h.in with current values
configure_file(${PROJECT_SOURCE_DIR}/source/config/Doxyfile.in ${PROJECT_SOURCE_DIR}/doc/doxygen/Doxyfile)

##TODO: when building documentation, set in environment using svnversion before calling doxygen
#if (NOT PROG_Svnversion-NOTFOUND)
#
#	# define some filenames
#	set( svnversion_file ${PROJECT_SOURCE_DIR}/include/openms_svnversion.txt) # file with just the revision number
#	set( svnversion_file_h ${PROJECT_SOURCE_DIR}/include/openms_svnversion.h) # file with #define directive, e.g. #define OPENMS_SVNREV 1234:2345M
	
	#set_source_files_properties(${PROJECT_SOURCE_DIR}/source/CONCEPT/VersionInfo.C)
	
	## ensure that the file is always rebuild
#  add_custom_target(svnversion_file_target ALL  ${PROG_Svnversion} -n ${PROJECT_SOURCE_DIR} > svnversion_file
#  									COMMAND "file(READ ${svnversion_file} ENV_OPENMS_REVISION)" ## read content of file into variable
#  									COMMAND "set( ENV{OPENMS_REVISION} ${ENV_OPENMS_REVISION} )" ## set env variable
#                   	COMMENT "Call svnversion in project directory to generate openms_svnversion.txt")
#	  
#  
#endif()

########################################################
###    BUILD the lib																 ###
########################################################

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

MACRO (QT4_WRAP_UI_OWN outfiles )
  QT4_EXTRACT_OPTIONS(ui_files ui_options ${ARGN})

  FOREACH (it ${ui_files}) 
    GET_FILENAME_COMPONENT(outfile ${it} NAME_WE)
    GET_FILENAME_COMPONENT(infile ${it} ABSOLUTE)
    SET(outfile ${PROJECT_SOURCE_DIR}/${directory}/ui_${outfile}.h)
    ADD_CUSTOM_COMMAND(OUTPUT ${outfile}
      COMMAND ${QT_UIC_EXECUTABLE}
      ARGS ${ui_options} -o ${outfile} ${infile}
      MAIN_DEPENDENCY ${infile})
    SET(${outfiles} ${${outfiles}} ${outfile})
  ENDFOREACH (it)

ENDMACRO (QT4_WRAP_UI_OWN)

set(OpenMS_sources  CACHE INTERNAL "This variable should hold all OpenMS sources at the end of the config step" )
#include(source/TEST/sources.cmake)                                   
include(source/KERNEL/sources.cmake)                                 
include(source/METADATA/sources.cmake)                               
include(source/ANALYSIS/QUANTITATION/sources.cmake)                  
include(source/ANALYSIS/SVM/sources.cmake)                           
include(source/ANALYSIS/MAPMATCHING/sources.cmake)
include(source/ANALYSIS/DECHARGING/sources.cmake)
include(source/ANALYSIS/ID/sources.cmake)
include(source/FILTERING/BASELINE/sources.cmake)
include(source/FILTERING/TRANSFORMERS/sources.cmake)
include(source/FILTERING/DATAREDUCTION/sources.cmake)
include(source/FILTERING/CALIBRATION/sources.cmake)
include(source/FILTERING/SMOOTHING/sources.cmake)
include(source/FILTERING/NOISEESTIMATION/sources.cmake)
include(source/FILTERING/ID/sources.cmake)
include(source/SYSTEM/sources.cmake)
include(source/APPLICATIONS/TOPP/sources.cmake)
include(source/APPLICATIONS/sources.cmake)
include(source/FORMAT/DB/sources.cmake)
include(source/FORMAT/HANDLERS/sources.cmake)
include(source/FORMAT/VALIDATORS/sources.cmake)
include(source/FORMAT/sources.cmake)
include(source/TRANSFORMATIONS/FEATUREFINDER/sources.cmake)
include(source/TRANSFORMATIONS/RAW2PEAK/sources.cmake)
include(source/MATH/STATISTICS/sources.cmake)
include(source/MATH/MISC/sources.cmake)
include(source/CONCEPT/sources.cmake)
include(source/VISUAL/DIALOGS/sources.cmake)
include(source/VISUAL/VISUALIZER/sources.cmake)
include(source/VISUAL/sources.cmake)
include(source/CHEMISTRY/sources.cmake)
include(source/COMPARISON/CLUSTERING/sources.cmake)
include(source/COMPARISON/SPECTRA/sources.cmake)
include(source/DATASTRUCTURES/sources.cmake)
#include(source/EXAMPLES/sources.cmake)
include(include/OpenMS/VISUAL/DIALOGS/UIC/sources.cmake)
include(include/OpenMS/VISUAL/UIC/sources.cmake)
include(include/OpenMS/VISUAL/DIALOGS/sources.cmake)
include(include/OpenMS/VISUAL/VISUALIZER/sources.cmake)
include(include/OpenMS/VISUAL/sources.cmake)
include(include/OpenMS/SYSTEM/sources.cmake)
include(include/OpenMS/APPLICATIONS/sources.cmake)
include(include/OpenMS/FORMAT/sources.cmake)


if (MINGW)
	ADD_LIBRARY(OpenMS STATIC ${OpenMS_sources} )
else()
	ADD_LIBRARY(OpenMS SHARED ${OpenMS_sources} )
endif()


# TODO track why the duplicate warnings are thrown for all (!) MOC sources
# Macro problem?
list(REMOVE_DUPLICATES OpenMS_sources)

ADD_LIBRARY(OpenMS SHARED ${OpenMS_sources} )

####### TOPP #########
set(TOPP_executables)
include(source/APPLICATIONS/TOPP/executables.cmake)

foreach(i ${TOPP_executables}) 
	add_executable(${i} source/APPLICATIONS/TOPP/${i}.C)
	target_link_libraries(${i} ${QT_LIBRARIES} ${CGAL_LIBRARY} ${GSL_LIBRARY} ${XERCESC_LIBRARY} /share/opt/x86_64_sl4/OpenMS_contrib-current-gcc_4.1.1/lib/libgslcblas.so /share/opt/x86_64_sl4/OpenMS_contrib-current-gcc_4.1.1/lib/svm.o OpenMS)
endforeach(i)
add_custom_target(TOPP)
add_dependencies(TOPP ${TOPP_executables})


##### UTILS ########
set(UTILS_executables)
include(source/APPLICATIONS/UTILS/executables.cmake)

foreach(i ${UTILS_executables})
	add_executable(${i} source/APPLICATIONS/UTILS/${i}.C)
	target_link_libraries(${i} ${QT_LIBRARIES} ${CGAL_LIBRARY} ${GSL_LIBRARY} ${XERCESC_LIBRARY} /share/opt/x86_64_sl4/OpenMS_contrib-current-gcc_4.1.1/lib/libgslcblas.so /share/opt/x86_64_sl4/OpenMS_contrib-current-gcc_4.1.1/lib/svm.o OpenMS)
endforeach(i)
add_custom_target(UTILS)
add_dependencies(UTILS ${UTILS_executables})


###### Examples ######
set(EXAMPLES_executables)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY source/EXAMPLES)
include(source/EXAMPLES/executables.cmake)
foreach(i ${EXAMPLES_executables})
  add_executable(${i} source/EXAMPLES/${i}.C)
  target_link_libraries(${i} ${QT_LIBRARIES} ${CGAL_LIBRARY} ${GSL_LIBRARY} ${XERCESC_LIBRARY} /share/opt/x86_64_sl4/OpenMS_contrib-current-gcc_4.1.1/lib/libgslcblas.so /share/opt/x86_64_sl4/OpenMS_contrib-current-gcc_4.1.1/lib/svm.o OpenMS)
endforeach(i)
add_custom_target(EXAMPLES)
add_dependencies(EXAMPLES ${EXAMPLES_executables})


##### Test Executables ####
set(TEST_executables)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY source/TEST)
include(source/TEST/executables.cmake)
foreach(i ${TEST_executables})
  add_executable(${i} source/TEST/${i}.C)
  target_link_libraries(${i} ${QT_LIBRARIES} ${CGAL_LIBRARY} ${GSL_LIBRARY} ${XERCESC_LIBRARY} /share/opt/x86_64_sl4/OpenMS_contrib-current-gcc_4.1.1/lib/libgslcblas.so /share/opt/x86_64_sl4/OpenMS_contrib-current-gcc_4.1.1/lib/svm.o OpenMS)
endforeach(i)
add_custom_target(TEST)
add_dependencies(TEST ${TEST_executables})


#### Tests #####
enable_testing()
add_subdirectory(source/TEST)

##### TOPPtests #####
add_subdirectory(source/TEST/TOPP)

                    
                    
########################################################
###    install/copy??																 ###
########################################################

# TODO: install (copy) directories?!
# ...


