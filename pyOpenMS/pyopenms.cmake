find_package(PythonInterp REQUIRED)

# find out python version info
execute_process(
     COMMAND
     ${PYTHON_EXECUTABLE} -c "import sys; print '%s.%s' % sys.version_info[:2]"
     OUTPUT_VARIABLE PY_VERSION
     OUTPUT_STRIP_TRAILING_WHITESPACE
)

MESSAGE(STATUS "found python ${PY_VERSION}")

# windows support restrecited to pyhton2.7 at the moment !
IF (WIN32)
    IF (NOT PY_VERSION STREQUAL "2.7")
        MESSAGE(STATUS "need python 2.7 on windows")
        RETURN()
    ENDIF()

    IF (NOT MSVC90)
        MESSAGE(STATUS "need visual c++ 2008 compiler for building python 2.7 extensions")
        RETURN()
    ENDIF()


    include(InstallRequiredSystemLibraries)
    SET(MSVCR90DLL ${MSVC90_CRT_DIR}/msvcr90.dll)
    IF (NOT EXISTS ${MSVCR90DLL})
        MESSAGE(STATUS "need visual c++ 2008 runtime (called vcredist)")
        RETURN()
    ENDIF()

ENDIF(WIN32)


MESSAGE(STATUS "Looking for cython")
find_program( CYTHON_EXECUTABLE NAMES cython )

SET(CYTHON-MISSING FALSE)
IF (DEFINED CYTHON_EXECUTABLE-NOTFOUND)
	SET(CYTHON-MISSING TRUE)
ENDIF()

IF (CYTHON-MISSING)
	MESSAGE(STATUS "Looking for cython - not found")
ELSE()
	MESSAGE(STATUS "Looking for cython - found")
ENDIF()

MESSAGE(STATUS "Looking for numpy")
execute_process(
     COMMAND
     ${PYTHON_EXECUTABLE} -c "import numpy"
     RESULT_VARIABLE NUMPY_MISSING
     ERROR_QUIET
     OUTPUT_QUIET
)

SET(NUMPY-MISSING TRUE)
IF( NUMPY_MISSING EQUAL 0)
  SET(NUMPY-MISSING FALSE)
ENDIF()
IF(NUMPY_MISSING)
	MESSAGE(STATUS "Looking for numpy - not found")
ELSE()
	MESSAGE(STATUS "Looking for numpy - found")
ENDIF()


IF (NUMPY-MISSING OR CYTHON-MISSING)
   RETURN()
ENDIF()


# copy files

FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/pyOpenMS)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/unittests)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/pyOpenMS/cython_code)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/pyOpenMS/cython_code/pxd)

FILE(GLOB _python_files "pyOpenMS/pyOpenMS/*.py")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyOpenMS)

FILE(GLOB _python_files "pyOpenMS/unittests/*")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/unittests)

FILE(GLOB _files "pyOpenMS/pyOpenMS/cython_code/*.py")
FILE(COPY ${_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyOpenMS/cython_code)

FILE(GLOB _python_files "pyOpenMS/pyOpenMS/cython_code/*.pyx")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyOpenMS/cython_code)

FILE(GLOB _python_files "pyOpenMS/pyOpenMS/cython_code/pxd/*.p*")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyOpenMS/cython_code/pxd)

FILE(COPY pyOpenMS/setup.py DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS)
FILE(COPY pyOpenMS/version.py DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS)
FILE(COPY pyOpenMS/run_nose.py DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS)
FILE(COPY pyOpenMS/build_zip_for_install.py DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS)
FILE(COPY pyOpenMS/main_for_installer.py DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS)

IF (WIN32)
    FILE(COPY ${MSVCR90DLL} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyOpenMS})
    SET(FOUND_XERCES FALSE)
    FOREACH(CONTRIB_PATH ${CONTRIB_DIR})
        IF (EXISTS ${CONTRIB_PATH}/lib/xerces-c_3_0.dll)
            FILE(COPY ${CONTRIB_PATH}/lib/xerces-c_3_0.dll DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyOpenMS)
            SET(FOUND_XERCES TRUE)
        ENDIF()
    ENDFOREACH()
    IF (NOT FOUND_XERCES)
        MESSAGE(STATUS "cound not find xerces dll in contrib dir")
        RETURN()
    ENDIF()
ENDIF()


# generate cython wrapper

MESSAGE(STATUS "Generate cython source file")
EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} build_cython_file.py 
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/pyOpenMS/cython_code
  ERROR_VARIABLE PYOK OUTPUT_QUIET)
MESSAGE(STATUS "Generate cython source file - done")


# run cython to generate c++ file
 
MESSAGE(STATUS "run cython to generate c++ file")
EXECUTE_PROCESS(COMMAND ${CYTHON_EXECUTABLE} -X boundscheck=False -X wraparound=False --cplus -o ../pyOpenMS.cpp pyOpenMS.pyx
WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/pyOpenMS/cython_code
ERROR_VARIABLE CYOK OUTPUT_QUIET)
MESSAGE(STATUS "run cython to generate c++ file - done")

# write variables for setup.py as python script

set(ENVPATH ${CMAKE_BINARY_DIR}/pyOpenMS/env.py)

FILE(WRITE ${ENVPATH} OPEN_MS_SRC="${CMAKE_SOURCE_DIR}" "\n" )
FILE(APPEND ${ENVPATH} OPEN_MS_BUILD_DIR="${CMAKE_BINARY_DIR}" "\n")

FILE(APPEND ${ENVPATH} OPEN_MS_CONTRIB_BUILD_DIRS=\")
FOREACH(CONTRIB_PATH ${CONTRIB_DIR})
	FILE(APPEND ${ENVPATH} ${CONTRIB_PATH} ";")
ENDFOREACH()
FILE(APPEND ${ENVPATH} "\"\n")

FILE(APPEND ${ENVPATH} QT_HEADERS_DIR="${QT_HEADERS_DIR}" "\n")
FILE(APPEND ${ENVPATH} QT_LIBRARY_DIR="${QT_LIBRARY_DIR}" "\n")
FILE(APPEND ${ENVPATH} QT_QTCORE_INCLUDE_DIR="${QT_QTCORE_INCLUDE_DIR}" "\n")
FILE(APPEND ${ENVPATH} MSVCR90DLL="${QT_QTCORE_INCLUDE_DIR}" "\n")
IF (WIN32)
    FILE(APPEND ${ENVPATH} OPEN_MS_LIB="${CMAKE_BINARY_DIR}/bin/OpenMS.dll" "\n")
ENDIF()

# create targets in makefile 

add_custom_target(pyopenms 
	COMMAND ${PYTHON_EXECUTABLE} setup.py build_ext --inplace
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
add_dependencies(pyopenms OpenMS)

add_custom_target(bdist_pyopenms 
	COMMAND ${PYTHON_EXECUTABLE} setup.py bdist  --format=zip
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
add_dependencies(bdist_pyopenms OpenMS)

add_custom_target(rpm_pyopenms 
	COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_rpm  
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
add_dependencies(rpm_pyopenms OpenMS)

add_custom_target(install_pyopenms 
	COMMAND ${PYTHON_EXECUTABLE} setup.py install 
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
add_dependencies(install_pyopenms OpenMS)

add_custom_target(pyopenms_installer 
	COMMAND ${PYTHON_EXECUTABLE} build_zip_for_install.py
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
add_dependencies(pyopenms_installer bdist_pyopenms)



    MESSAGE(STATUS )
    MESSAGE(STATUS ---------------------------------------------------------------- )
    MESSAGE(STATUS )
    MESSAGE(STATUS "Additional build targets:" )
    MESSAGE(STATUS )
    MESSAGE(STATUS "   pyopenms           -- builds python extension module")
    MESSAGE(STATUS "   bdist_pyopenms     -- builds binary distribution package")
IF (NOT WIN32)
    MESSAGE(STATUS "   rpm_pyopenms       -- builds rpm distribution package")
ENDIF()
    MESSAGE(STATUS "   install_pyopenms   -- global install of pyopenms ")
    MESSAGE(STATUS "   pyopenms_installer -- builds installer as python zip package")
    MESSAGE(STATUS )
    MESSAGE(STATUS ---------------------------------------------------------------- )
    MESSAGE(STATUS )


enable_testing()
add_test(NAME test_pyopenms 
         COMMAND ${PYTHON_EXECUTABLE} run_nose.py
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS 
        )
IF(NOT WIN32)
    set_tests_properties(test_pyopenms PROPERTIES ENVIRONMENT 
                         "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib")
    set_target_properties(pyopenms PROPERTIES ENVIRONMENT 
                         "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib")
ENDIF()
