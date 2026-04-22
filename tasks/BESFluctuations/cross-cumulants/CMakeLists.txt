include_directories ("${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_BINARY_DIR}/include")

add_executable (BES-cross-cumulants BES-cross-cumulants.cpp)
target_link_libraries (BES-cross-cumulants ThermalFIST)
set_property(TARGET BES-cross-cumulants PROPERTY FOLDER "fist-sampler")