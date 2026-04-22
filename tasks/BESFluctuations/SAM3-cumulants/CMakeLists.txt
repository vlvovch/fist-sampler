include_directories ("${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_BINARY_DIR}/include")

add_executable (BES-SAM3-cumulants BES-SAM3-cumulants.cpp)
target_link_libraries (BES-SAM3-cumulants ThermalFIST)
set_property(TARGET BES-SAM3-cumulants PROPERTY FOLDER "fist-sampler")
