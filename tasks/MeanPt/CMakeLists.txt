include_directories ("${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_BINARY_DIR}/include")

add_executable (MeanPtCalculation ../MeanPt/MeanPtCalculation.cpp)
#target_link_libraries (CooperFryeSampler ThermalFIST ThermalFISTExtended sample_statistics)
target_link_libraries (MeanPtCalculation ThermalFIST)
set_property(TARGET MeanPtCalculation PROPERTY FOLDER "fist-sampler")
