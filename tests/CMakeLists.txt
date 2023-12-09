cmake_minimum_required(VERSION 3.5)
project(tests VERSION 0.1.0)

#add_definitions(-DTEST=42)

CPMAddPackage(
        NAME Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.3.2
        EXCLUDE_FROM_ALL YES
        GIT_SHALLOW YES
)

include(CTest)
enable_testing()

include_directories(${tests_SOURCE_DIR}/include)

message(${SOURCES})

add_library(kstest STATIC ksTest.cpp)
target_include_directories(kstest PUBLIC ${tests_SOURCE_DIR}/include)
link_libraries(Catch2::Catch2WithMain dpm_validation kstest)

add_subdirectory(cpu)

if (GPU)
    add_subdirectory(gpu)
endif ()
