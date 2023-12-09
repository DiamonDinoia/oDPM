link_libraries(ODPMLIB)

if (GPU)
    file(GLOB_RECURSE TEST_SOURCES "./*.cpp")
    set_source_files_properties(${TEST_SOURCES} PROPERTIES LANGUAGE CUDA)
endif ()

add_executable(validation_tests validation_tests.cpp)

add_executable(geometry_tests geometry_tests.cpp)

add_executable(electron_tracking_tests electron_tracking_tests.cpp)

add_executable(photon_simulation_tests photon_simulation_tests.cpp)

add_executable(simulation_tests simulation_tests.cpp)

add_executable(table_tests table_tests.cpp)

add_test(NAME interactions COMMAND validation_tests)
add_test(NAME geometry COMMAND geometry_tests)
add_test(NAME electron-tracking COMMAND electron_tracking_tests)
add_test(NAME photon-simulation COMMAND photon_simulation_tests)
add_test(NAME simulation COMMAND simulation_tests)
add_test(NAME tables COMMAND table_tests)
