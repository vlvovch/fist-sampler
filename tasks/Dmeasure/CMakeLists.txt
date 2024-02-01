include_directories ("${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_BINARY_DIR}/include")

add_executable (DmeasureCalculationLHC DmeasureCalculationLHC.cpp)
#target_link_libraries (CooperFryeSampler ThermalFIST ThermalFISTExtended sample_statistics)
target_link_libraries (DmeasureCalculationLHC ThermalFIST)
set_property(TARGET DmeasureCalculationLHC PROPERTY FOLDER "fist-sampler")
