# CMake installation file
# Created 2008/10 during OpenMS retreat by CB && AB

PROJECT("OpenMS")

########################################################
###    manual entries (edit this for new release)    ###
########################################################

set(CF_PACKAGE_VERSION "1.3")

########################################################
###    end manual																	   ###
########################################################

# Heart of the BUILD system : only edit when you know what you´re doing (we don´t)
# quick manual for most commands: http://www.cmake.org/cmake/help/cmake2.6docs.html

cmake_minimum_required(VERSION 2.6)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

## define some directories
if (CMAKE_INSTALL_PREFIX)
	set(DIR_SHARE ${CMAKE_INSTALL_PREFIX}/share_${CF_PACKAGE_VERSION}/ CACHE INTERNAL)
else (CMAKE_INSTALL_PREFIX)
	set(DIR_SHARE ${PROJECT_SOURCE_DIR}/share/ CACHE INTERNAL)
endif (CMAKE_INSTALL_PREFIX)


AUX_SOURCE_DIRECTORY(source/ANALYSIS ANALYSIS_SOURCES)
AUX_SOURCE_DIRECTORY(source/APPLICATIONS APPLICATIONS_SOURCES)
AUX_SOURCE_DIRECTORY(source/CHEMISTRY CHEMISTRY_SOURCES)
AUX_SOURCE_DIRECTORY(source/COMPARISON COMPARISON_SOURCES)
AUX_SOURCE_DIRECTORY(source/CONCEPT CONCEPT_SOURCES)
AUX_SOURCE_DIRECTORY(source/DATASTRUCTURES DATASTRUCTURES_SOURCES)
AUX_SOURCE_DIRECTORY(source/FILTERING FILTERING_SOURCES)
AUX_SOURCE_DIRECTORY(source/FORMAT FORMAT_SOURCES)
AUX_SOURCE_DIRECTORY(source/KERNEL KERNEL_SOURCES)
AUX_SOURCE_DIRECTORY(source/MATH MATH_SOURCES)
AUX_SOURCE_DIRECTORY(source/METADATA METADATA_SOURCES)
AUX_SOURCE_DIRECTORY(source/SYSTEM SYSTEM_SOURCES)
AUX_SOURCE_DIRECTORY(source/TRANSFORMATIONS TRANSFORMATIONS_SOURCES)
AUX_SOURCE_DIRECTORY(source/VISUAL/DIALOGS VISUAL_SOURCES)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include/)

########################################################
###    compiler flags																 ###
########################################################

## general
#TODO confirm that these flags are general
add_definitions("-Wall -ansi -pedantic -Wno-long-long -Wno-non-virtual-dtor")

## platform dependent compiler flags:
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-fPIC" WITH_FPIC)
if (WITH_FPIC)
	add_definitions(-fPIC)
endif()



########################################################
###    external libs																 ###
########################################################


FIND_PACKAGE(Qt4 REQUIRED)
SET(QT_USE_QTXML 1)

# bla add others, only QtCore and QtGui are enabled by default
INCLUDE(${QT_USE_FILE})

INCLUDE(UseQt4)



#### bla
MACRO (QT4_WRAP_UI_OWN outfiles dest_dir)
    QT4_EXTRACT_OPTIONS(ui_files ui_options ${ARGN})
    FOREACH (it ${ui_files}) 
      GET_FILENAME_COMPONENT(outfile ${it} NAME_WE)
      GET_FILENAME_COMPONENT(infile ${it} ABSOLUTE)
      SET(outfile ${dest_dir}/${outfile}Template.h)
			#message(${QT_UIC_EXECUTABLE})
			message("Generating uic header of file: " ${outfile})
			execute_process(COMMAND ${QT_UIC_EXECUTABLE} ${ui_options} -o ${outfile} ${infile})
      SET(${outfiles} ${${outfiles}} ${outfile})
    ENDFOREACH (it)
ENDMACRO (QT4_WRAP_UI_OWN)

ADD_SUBDIRECTORY(source/VISUAL)
ADD_SUBDIRECTORY(source/VISUAL/DIALOGS)

# TODO check for Cmake custom flag for this 
set(CONTRIB_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/contrib/)
set(CONTRIB_LIB_DIR ${PROJECT_SOURCE_DIR}/contrib/)

# TODO check whether this is necessary
FIND_PATH(GSL_INCLUDE_DIR
    NAMES gsl/gsl_cdf.h gsl/gsl_randist.h
    PATHS ${CONTRIB_INCLUDE_DIR}
    PATH_SUFFIXES include
    DOC "GSL header include dir"
    )
MESSAGE(STATUS "Found gsl headers at: " ${GSL_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIR})

# TODO check whether this is necessary
FIND_LIBRARY(GSL_GSL_LIBRARY
    NAMES gsl libgsl
    PATHS  ${CONTRIB_LIB_DIR}
    PATH_SUFFIXES lib
    DOC "GSL library dir" 
		)
MESSAGE(STATUS "Found gsl library at: " ${GSL_GSL_LIBRARY})

ADD_LIBRARY(OpenMS SHARED 
						${ANALYSIS_SOURCES} 
						${APPLICATIONS_SOURCES}
						{CHEMISTRY_SOURCES}
						${COMPARISON_SOURCES}
						${CONCEPT_SOURCES}
						${DATASTRUCTURES_SOURCES}
						${FILTERING_SOURCES}
						${FORMAT_SOURCES}
						${KERNEL_SOURCES}
						${MATH_SOURCES}
						${METADATA_SOURCES}
						${SYSTEM_SOURCES}
						${TRANSFORMATIONS_SOURCES}
						${VISUAL_DIALOGS_UI_SOURCES}
						${VISUAL_SOURCES}
						)


########################################################
###    useful programms															 ###
########################################################

find_program(PROG_Svnversion             
             svnversion
             PATHS "c:\programme\subversion\bin\" "c:\program files\subversion\bin"   ### additional search paths (along with $PATH)
             DOC "svnversion executable which helps in determining the svn revision when building TOPP tools"

if (PROG_Svnversion-NOTFOUND)
	message(STATUS "The programm svnversion could not be found. SVN-revision information will not be available! Add the location of svnversion(.exe) to your PATH environment variable if you require SVN-revision.")
endif()

########################################################
###    configure config.h														 ###
########################################################


# user DEBUG flag (needs a cmake flag/cache entry, so we define it empty which causes it to show up in the wizard)
set(CF_OPENMS_DEBUG)

# can be set automatically
set(CF_OPENMS_DATA_PATH ${DIR_SHARE}/OpenMS/ CACHE INTERNAL)

# check for Microsoft Visual Studio compiler
if (MSVC)
	set(CF_OPENMS_COMPILER_MSVC "1" CACHE INTERNAL);
endif()
# check for G++
if (CMAKE_COMPILER_IS_GNUCXX)
	set(CF_OPENMS_COMPILER_GXX "1" CACHE INTERNAL);
endif()

INCLUDE(TestBigEndian)
TEST_BIG_ENDIAN(CF_OPENMS_BIG_ENDIAN)

# check 32/64 bit architecture
if (CMAKE_SIZEOF_VOID_P MATCHES "8")
	set(CF_OPENMS_64BIT_ARCHITECTURE, "1")
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

CHECK_INCLUDE_FILE_CXX("unistd.h" CF_OPENMS_HAS_UNISTD_H)
CHECK_INCLUDE_FILE_CXX("process.h" CF_OPENMS_HAS_PROCESS_H)

CHECK_INCLUDE_FILE_CXX("time.h" CF_OPENMS_HAS_TIME_H)
CHECK_INCLUDE_FILE_CXX("sys/types.h" CF_OPENMS_HAS_SYS_TYPES_H)
CHECK_INCLUDE_FILE_CXX("sys/times.h" CF_OPENMS_HAS_SYS_TIMES_H)
CHECK_INCLUDE_FILE_CXX("sys/time.h"  CF_OPENMS_HAS_SYS_TIME_H)

# user flag (needs a cmake flag/cache entry, so we define it empty which causes it to show up in the wizard)
set(CF_OPENMS_ENABLE_EXTERNALMEMORY)

include(CheckFunctionExists)
## in MinGW we have the signal.h header, but no kill() as in Linux, so we need to check for the kill() function 
CHECK_FUNCTION_EXISTS("kill" CF_OPENMS_HAS_KILL)
CHECK_FUNCTION_EXISTS("sysconf" CF_OPENMS_HAS_SYSCONF)

# user flag with default "QMYSQL" (needs a cmake flag/cache entry, so we define it empty which causes it to show up in the wizard)
set(CF_DB_PLUGIN "QMYSQL")


## replace any variables in config.h.in with current values
configure_file(${PROJECT_SOURCE_DIR}/source/config/config.h.in ${PROJECT_SOURCE_DIR}/include/config.h)

########################################################
###    configure Doxyfile														 ###
########################################################

## CF_PACKAGE_VERSION (already defined above in manual section)

set(CF_OPENMS_PATH ${PROJECT_BINARY_DIR})

## replace any variables in config.h.in with current values
configure_file(${PROJECT_SOURCE_DIR}/source/config/Doxyfile.in ${PROJECT_SOURCE_DIR}/doc/doxygen/Doxyfile)

##TODO: when building documentation, set in environment using svnversion before calling doxygen
if (NOT PROG_Svnversion-NOTFOUND)

	# define some filenames
	set( svnversion_file ${PROJECT_SOURCE_DIR}/include/openms_svnversion.txt) # file with just the revision number
	set( svnversion_file_h ${PROJECT_SOURCE_DIR}/include/openms_svnversion.h) # file with #define directive, e.g. #define OPENMS_SVNREV 1234:2345M
	
	#set_source_files_properties(${PROJECT_SOURCE_DIR}/source/CONCEPT/VersionInfo.C)
	
	## ensure that the file is always rebuild
  add_custom_target(svnversion_file_target ALL  ${PROG_Svnversion} -n ${PROJECT_SOURCE_DIR} > svnversion_file
                   	COMMENT "Call svnversion in project directory to generate openms_svnversion.txt")
	 
 
  file(READ ${svnversion_file} ENV_OPENMS_REVISION) ## read content of file into variable
  set( ENV{OPENMS_REVISION} ${ENV_OPENMS_REVISION} )## set env variable
endif()

                    
                    
########################################################
###    install/copy??																 ###
########################################################

# TODO: install (copy) directories?!
# ...
