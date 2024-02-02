project(FISTSampler)

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

set (CMAKE_CXX_STANDARD 11)

# The version number.
set (FistSampler_VERSION_MAJOR 1)
set (FistSampler_VERSION_MINOR 0)

#include_directories ("${PROJECT_SOURCE_DIR}/include" "${PROJECT_BINARY_DIR}/include")

if(NOT TARGET ThermalFIST)
  add_subdirectory(Thermal-FIST)
endif()

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/FistSamplerConfig.h.in"
  "${PROJECT_BINARY_DIR}/include/FistSamplerConfig.h"
  )

include_directories ("${PROJECT_SOURCE_DIR}/include" 
  "${PROJECT_BINARY_DIR}/include")

add_executable (FISTSampler FistSampler.cpp include/HypersurfaceReader.h include/FistSamplerParameters.h  include/FistSamplerHelperFunctions.h)
#target_link_libraries (CooperFryeSampler ThermalFIST ThermalFISTExtended sample_statistics)
target_link_libraries (FISTSampler ThermalFIST)
set_property(TARGET FISTSampler PROPERTY FOLDER "fist-sampler")

add_subdirectory(tasks/Dmeasure)
add_subdirectory(tasks/MeanPt)