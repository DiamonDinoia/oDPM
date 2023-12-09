link_libraries(GPUODPMLIB)

add_executable(gpu_table_tests table_tests.cu)
add_test(NAME gpu_tables COMMAND gpu_table_tests)

add_executable(gpu_geometry_tests geometry_tests.cu)
add_test(NAME gpu_geometry COMMAND geometry_tests)

add_executable(gpu_photon_simulations photon_simulations_tests.cu)
add_test(NAME gpu_photon_simulation COMMAND gpu_photon_simulations)

add_executable(gpu_simulation_tests simulation_tests.cu)
add_test(NAME gpu_simulation COMMAND gpu_simulation_tests)

add_executable(gpu_electron_tracking electron_tracking_tests.cu)
add_test(NAME gpu-electron-tracking COMMAND gpu_electron_tracking)
