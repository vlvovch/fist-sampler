include_directories ("${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_BINARY_DIR}/include")

add_executable (BES-alphas-KpLaQ BES-alphas-KpLaQ.cpp)
#target_link_libraries (CooperFryeSampler ThermalFIST ThermalFISTExtended sample_statistics)
target_link_libraries (BES-alphas-KpLaQ ThermalFIST)
set_property(TARGET BES-alphas-KpLaQ PROPERTY FOLDER "fist-sampler")