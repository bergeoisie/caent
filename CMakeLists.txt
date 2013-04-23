cmake_minimum_required (VERSION 2.6)
project (caent)
 
find_package(armadillo)
find_package(boost)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/caentConfig.h.in"
  "${PROJECT_BINARY_DIR}/caentConfig.h"
  )
 
# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")

set(CMAKE_CXX_FLAGS "-g -Wall -O2")

add_library(caenthelper caenthelper.cpp)

# add the executable
add_executable(probrange probrange.cpp)

target_link_libraries(probrange armadillo caenthelper)