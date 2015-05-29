
# @author Shin'ichiro Nakaoka

include(CheckFunctionExists)

if(UNIX)
    if(NOT LAPACK_LIBRARY_DIRS)

      find_library(
        LAPACK_LIBRARY lapack
        PATHS /usr/lib/atlas /usr/lib/atlas-base/atlas
        NO_DEFAULT_PATH)

      if(NOT LAPACK_LIBRARY)
        find_library(LAPACK_LIBRARY lapack)
      endif()

      if(LAPACK_LIBRARY)
        get_filename_component(LAPACK_LIBRARY_DIRS ${LAPACK_LIBRARY} PATH)
        message(STATUS "detected ${LAPACK_LIBRARY}")
      endif()

    endif()

    if(NOT LAPACK_LIBRARIES AND LAPACK_LIBRARY_DIRS)

      find_library(
        BLAS_LIBRARY cblas
        PATHS ${LAPACK_LIBRARY_DIRS}
        NO_DEFAULT_PATH)

      if(NOT BLAS_LIBRARY)
        find_library(
          BLAS_LIBRARY blas
          PATHS ${LAPACK_LIBRARY_DIRS}
          NO_DEFAULT_PATH)
      endif()

      find_library(
        LAPACK_LIBRARY lapack
        PATH ${LAPACK_LIBRARY_DIRS}
        NO_DEFAULT_PATH)

      if(BLAS_LIBRARY AND LAPACK_LIBRARY)
        list(APPEND LAPACK_LIBRARIES ${BLAS_LIBRARY} ${LAPACK_LIBRARY})
      endif()

      find_library(
        G2C_LIBRARY g2c
        PATH ${LAPACK_LIBRARY_DIRS}
        NO_DEFAULT_PATH)

      if(G2C_LIBRARY)
        list(APPEND LAPACK_LIBRARIES ${G2C_LIBRARY})
      endif()

    endif()
endif(UNIX)

if(WIN32)
  set(LAPACK_LIBRARIES optimized lapack debug lapackd optimized blas debug blasd optimized libf2c debug libf2cd)
  if( ${CMAKE_GENERATOR} MATCHES "Win64$" )
	if(NOT LAPACK_TOP_DIR)
		  find_path(
		  LAPACK_TOP_DIR 
		  NAMES LIB/x64/clapack.lib
		  PATHS $ENV{HOMEDRIVE}/Program Files $ENV{HOMEDRIVE}/
		  PATH_SUFFIXES CLAPACK CLAPACK3.1.1 CLAPACK-3.1.1
		  DOC "the top directory of clapack")
	  endif()
	  if(LAPACK_TOP_DIR)
		set(LAPACK_LIBRARY_DIRS ${LAPACK_TOP_DIR}/LIB/x64)
	  endif()
  else()
	  if(NOT LAPACK_TOP_DIR)
		  find_path(
		  LAPACK_TOP_DIR 
		  NAMES LIB/Win32/clapack.lib
		  PATHS $ENV{HOMEDRIVE}/Program Files $ENV{HOMEDRIVE}/
		  PATH_SUFFIXES CLAPACK CLAPACK3.1.1 CLAPACK-3.1.1
		  DOC "the top directory of clapack")
	  endif()
	  if(LAPACK_TOP_DIR)
		set(LAPACK_LIBRARY_DIRS ${LAPACK_TOP_DIR}/LIB/Win32)
	  endif()
  endif()
endif(WIN32)
  

if(NOT LAPACK_INCLUDE_DIRS)
if(UNIX)
  if(BLAS_LIBRARY)
    get_filename_component(BLAS_LIBRARY_DIR ${BLAS_LIBRARY} PATH)
    string(REGEX REPLACE "/lib/?" "" BLAS_DIR ${BLAS_LIBRARY_DIR})
    find_file(CBLAS_H_FILE cblas.h
      PATHS ${BLAS_DIR}/include /usr/local/include /usr/include)
    if(CBLAS_H_FILE)
      get_filename_component(BLAS_INCLUDE_DIR ${CBLAS_H_FILE} PATH)
      list(APPEND LAPACK_INCLUDE_DIRS ${BLAS_INCLUDE_DIR})
    endif()
  endif()

  if(LAPACK_LIBRARY)
    get_filename_component(LAPACK_LIBRARY_DIR ${LAPACK_LIBRARY} PATH)
    string(REGEX REPLACE "/lib/?" "" LAPACK_DIR ${LAPACK_LIBRARY_DIR})
    find_file(CLAPACK_H_FILE clapack.h
      PATHS ${LAPACK_DIR}/include /usr/local/include /usr/include)
    if(CLAPACK_H_FILE)
      get_filename_component(LAPACK_INCLUDE_DIR ${CLAPACK_H_FILE} PATH)
      list(APPEND LAPACK_INCLUDE_DIRS ${LAPACK_INCLUDE_DIR})
    endif()
  endif()
endif(UNIX)

if(WIN32)
  if(LAPACK_TOP_DIR)
    find_file(CLAPACK_H_FILE clapack.h
      PATHS ${LAPACK_TOP_DIR}/INCLUDE )
    mark_as_advanced(CLAPACK_H_FILE)
    if(CLAPACK_H_FILE)
      set( LAPACK_INCLUDE_DIRS ${LAPACK_TOP_DIR}/INCLUDE )
    endif()
  endif()    
endif(WIN32)

endif()

if(LAPACK_LIBRARIES)
  set(LAPACK_FOUND TRUE)

  set(CMAKE_REQUIRED_INCLUDES ${LAPACK_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES})

  CHECK_FUNCTION_EXISTS(clapack_dgesv clapack_dgesv_found)
  if(clapack_dgesv_found)
    set(USE_CLAPACK_INTERFACE TRUE)
    message(STATUS "use C interface for Lapack")
  endif()

  CHECK_FUNCTION_EXISTS(cblas_ddot cblas_ddot_found)
  if(cblas_ddot_found)
    set(USE_CBLAS_INTERFACE TRUE)
    message(STATUS "use c interaface for Blas")
  endif()
endif()


set(LAPACK_INCLUDE_DIRS ${LAPACK_INCLUDE_DIRS} CACHE PATH "Directories containing Lapack header files")
set(LAPACK_LIBRARY_DIRS ${LAPACK_LIBRARY_DIRS} CACHE PATH "Directories containing Lapack library files")
set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE PATH "Lapack library files to link")


if(LAPACK_FOUND)
  if(NOT Lapack_FIND_QUIETLY)
    message(STATUS "Found ${LAPACK_LIBRARIES} in ${LAPACK_LIBRARY_DIRS}")
  endif()
else()
  if(NOT Lapack_FIND_QUIETLY)
    if(Lapack_FIND_REQUIRED)
      message(FATAL_ERROR "Blas and Lapack required, please specify it's location.")
    endif()
  endif()
endif()
