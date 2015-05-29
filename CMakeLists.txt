cmake_minimum_required(VERSION 2.8)

project(USAC)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

option(USAC_BUILD_TEST "Build the USAC test program." OFF)

add_definitions(-DLIBCONFIGXX_STATIC)
find_package(Config++)
find_package(Lapack)

add_subdirectory(src)
if(USAC_BUILD_TEST)
	add_subdirectory(RunSingleTest)
endif(USAC_BUILD_TEST)
