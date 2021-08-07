project(CooperFryeSampler)

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

set (CMAKE_CXX_STANDARD 11)

#include_directories ("${PROJECT_SOURCE_DIR}/include" "${PROJECT_BINARY_DIR}/include")

if(NOT TARGET ThermalFIST)
  add_subdirectory(Thermal-FIST)
endif()

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/CooperFryeSamplerConfig.h.in"
  "${PROJECT_BINARY_DIR}/include/CooperFryeSamplerConfig.h"
  )

include_directories ("${PROJECT_SOURCE_DIR}/include" 
  "${PROJECT_BINARY_DIR}/include")

add_executable (CooperFryeSampler CooperFryeSampler.cpp include/HypersurfaceReader.h include/CooperFryeSamplerParameters.h  include/CooperFryeSamplerHelperFunctions.h)
#target_link_libraries (CooperFryeSampler ThermalFIST ThermalFISTExtended sample_statistics)
target_link_libraries (CooperFryeSampler ThermalFIST)
set_property(TARGET CooperFryeSampler PROPERTY FOLDER "research/CooperFryeSampler")
