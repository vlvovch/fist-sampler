include_directories ("${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_BINARY_DIR}/include")

add_executable (DmeasureCalculation DmeasureCalculation.cpp)
#target_link_libraries (CooperFryeSampler ThermalFIST ThermalFISTExtended sample_statistics)
target_link_libraries (DmeasureCalculation ThermalFIST)
set_property(TARGET DmeasureCalculation PROPERTY FOLDER "fist-sampler")
