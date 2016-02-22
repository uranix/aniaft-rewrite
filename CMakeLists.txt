cmake_minimum_required(VERSION 3.0)

project(aniaft)

set(aniaft_SOURCES
	boundary.cpp
	mesh.cpp
	tree.cpp
	tria.cpp
	)

add_library(aniaft SHARED ${aniaft_SOURCES})

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Werror -Wextra -g -funroll-loops")
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Werror -Wextra -g -funroll-loops")
endif()

find_package(SWIG 3.0 REQUIRED)
include(${SWIG_USE_FILE})

find_package(PythonLibs 2.7 REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

SET_SOURCE_FILES_PROPERTIES(aniaft.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(aniaft.i PROPERTIES SWIG_FLAGS "-Wall")
SWIG_ADD_MODULE(aniaft python aniaft.i)
SWIG_LINK_LIBRARIES(aniaft ${PYTHON_LIBRARIES} aniaft)
