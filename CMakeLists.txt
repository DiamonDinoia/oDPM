if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
    message(STATUS "AltMin requires an out-of-source build.")
    message(STATUS "Please remove these files from ${CMAKE_CURRENT_BINARY_DIR} first:")
    message(STATUS "CMakeCache.txt")
    message(STATUS "CMakeFiles")
    message(STATUS "Once these files are removed, create a separate directory")
    message(STATUS "and run CMake from there")
    message(FATAL_ERROR "in-source build detected")
endif ()

cmake_minimum_required(VERSION 3.19)
project(oDPM VERSION 0.1.0 LANGUAGES CXX C)

option(GENERATE_DATA "Generates required tables" OFF)
option(ENABLE_LTO "enables link time optimization" ON)
option(USE_AGGRESSIVE_OPTIMIZATIONS "enables bitwise arithmetic to reduce branches" ON)
option(ENABLE_VALIDATION "enables validation" ON)
option(PYTHON_BINDINGS "enables python bindings" ON)

if (SKBUILD)
    set(PYTHON_BINDINGS ON)
endif ()

if (PYTHON_BINDINGS)
    set(ENABLE_VALIDATION ON)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif ()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(DCMAKE_RULE_MESSAGES OFF)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(LIB_NAME oDPMlib)
set(EXE_NAME oDPM)

set(DATA_DIR ${PROJECT_SOURCE_DIR}/g4)
set(GEANT4_INSTALL_DATADIR ${PROJECT_SOURCE_DIR}/geant4)
file(MAKE_DIRECTORY ${DATA_DIR})
file(MAKE_DIRECTORY ${GEANT4_INSTALL_DATADIR})
add_definitions(-DDATA_DIR="${DATA_DIR}")

include(SystemPackagesConfig)

# USING CPM TO HANDLE DEPENDENCIES
set(CPM_DOWNLOAD_VERSION 0.38.5)
include(CpmConfig)
include(FlagsConfig)
add_definitions(-DABSL_ALLOCATOR_NOTHROW=1)

#
#
# -----------------------------------------------------------------------------
# Colored output
# -----------------------------------------------------------------------------

# From: https://medium.com/@alasher/colored-c-compiler-output-with-ninja-clang-gcc-10bfe7f2b949
option(FORCE_COLORED_OUTPUT "Always produce ANSI-colored output (GNU/Clang only)." FALSE)

if (${FORCE_COLORED_OUTPUT})
    if (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR (CMAKE_CXX_COMPILER MATCHES "g\\+\\+"))
        add_compile_options("-fdiagnostics-color=always")
    elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR (CMAKE_CXX_COMPILER MATCHES "clang\\+\\+"))
        add_compile_options("-fcolor-diagnostics")
    endif ()
endif ()

find_program(CCACHE_PROGRAM ccache)
if (CCACHE_PROGRAM)
    set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    set(CMAKE_CUDA_COMPILER_LAUNCHER "${CCACHE_PROGRAM}") # CMake 3.9+
endif ()


if (USE_AGGRESSIVE_OPTIMIZATIONS)
    add_definitions(-DUSE_AGGRESSIVE_OPTIMIZATIONS)
endif ()

include(LinkTimeOptimizationConfig)

include(DependenciesConfig)

set(DICOM_INCLUDE ${geant4_SOURCE_DIR}/examples/extended/medical/DICOM/include)
file(GLOB DICOM_SOURCES ${geant4_SOURCE_DIR}/examples/extended/medical/DICOM/src/*.cc)

set(DICOM2_INCLUDE ${geant4_SOURCE_DIR}/examples/extended/medical/DICOM2/include)
file(GLOB DICOM_SOURCES2 ${geant4_SOURCE_DIR}/examples/extended/medical/DICOM2/src/*.cc)

set(DICOM2_DATA_DIR ${geant4_SOURCE_DIR}/examples/extended/medical/DICOM2/data)
set(DICOM_DATA_DIR ${geant4_SOURCE_DIR}/examples/extended/medical/DICOM)

add_library(g4Dicom STATIC ${DICOM_SOURCES} ${DICOM_SOURCES2})
target_include_directories(g4Dicom PUBLIC ${DICOM_INCLUDE} ${DICOM2_INCLUDE})


file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dicomReader)
file(COPY ${DICOM_DATA_DIR}/dicomReader/include DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/dicomReader)
file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/dicomReader/include/DicomFileMgr.hh)
file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/dicomReader/include/DicomFileCT.hh)

set(DICOM_READER_INCLUDE ${CMAKE_CURRENT_BINARY_DIR}/dicomReader/include)

file(GLOB DICOM_READER_SOURCES ${DICOM_DATA_DIR}/dicomReader/src/*.cc)

file(GLOB DICOM_HEADERS ${DICOM_DATA_DIR}/dicomReader/include/*.hh)
message("DICOM_HEADERS: ${DICOM_HEADERS}")
message("DICOM_READER_INCLUDE: ${DICOM_READER_INCLUDE}")
message("DICOM_READER_SOURCES: ${DICOM_READER_SOURCES}")

set(CUSTOM_DICOM_READER_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/lib/dicomReader/include)
file(GLOB CUSTOM_DICOM_READER_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/lib/dicomReader/src/*.cc)
add_library(g4DicomReader STATIC ${DICOM_READER_SOURCES} ${CUSTOM_DICOM_READER_SOURCES})
target_include_directories(g4DicomReader PUBLIC ${DICOM_READER_INCLUDE} ${CUSTOM_DICOM_READER_INCLUDE})
target_include_directories(g4DicomReader PRIVATE
        ${dcmtk_SOURCE_DIR}/dcmdata/include
        ${dcmtk_SOURCE_DIR}/ofstd/include
        ${dcmtk_SOURCE_DIR}/oflog/include
        ${DCMTK_BINARY_DIR}/config/include/
)

target_link_libraries(g4DicomReader PUBLIC DCMTK config G4run-static G4zlib-static dcmdata G4geomtext-static)

target_link_libraries(g4Dicom PUBLIC g4DicomReader)
target_compile_definitions(g4Dicom PUBLIC -DG4_DCMTK)
list(APPEND LIBRARIES g4Dicom)

include(LanguageStandardConfig)

# finding all source files
file(GLOB_RECURSE SOURCES "${oDPM_SOURCE_DIR}/src/*.cpp")

# adding all the include directories
list(APPEND INCLUDE_DIRS ${oDPM_SOURCE_DIR}/include ${oDPM_SOURCE_DIR}/include/tables)

add_library(${LIB_NAME} OBJECT ${SOURCES})
target_link_libraries(${LIB_NAME} PUBLIC ${LIBRARIES})
target_include_directories(${LIB_NAME} PUBLIC ${INCLUDE_DIRS})
enable_lto(${LIB_NAME})
target_add_compilation_flags(${LIB_NAME})
set_target_properties(${LIB_NAME} PROPERTIES CXX_STANDARD ${ODPM_CXX_STANDARD})
target_compile_features(${LIB_NAME} PUBLIC cxx_std_${ODPM_CXX_STANDARD})

target_compile_definitions(${LIB_NAME} PUBLIC -DG4DATASET_DIR="${GEANT4_INSTALL_DATADIR}")
target_compile_definitions(${LIB_NAME} PUBLIC -DINPUT_DATA_DIR="${dpmg4cpp_SOURCE_DIR}/SBTables/data")
target_compile_definitions(${LIB_NAME} PUBLIC -DOUTPUT_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
target_compile_definitions(${LIB_NAME} PRIVATE -DDICOM2_DATA_DIR="${DICOM2_DATA_DIR}")
target_compile_definitions(${LIB_NAME} PRIVATE -DDICOM_DATA_DIR="${DICOM_DATA_DIR}")


if (ENABLE_VALIDATION)
    add_subdirectory(validation EXCLUDE_FROM_ALL)
    add_executable(validationExe validation.cpp)
    target_link_libraries(validationExe DPMValidation)
endif ()

function(print_all_targets DIR)
    get_property(TGTS DIRECTORY "${DIR}" PROPERTY BUILDSYSTEM_TARGETS)
    foreach(TGT IN LISTS TGTS)
        if (NOT(${TGT} IN_LIST LIBRARIES))
            set_target_properties(${TGT} PROPERTIES EXCLUDE_FROM_ALL TRUE EXCLUDE_FROM_DEFAULT_BUILD TRUE)
#            message(STATUS "Target: ${TGT}")
        endif()
    endforeach()
    get_property(SUBDIRS DIRECTORY "${DIR}" PROPERTY SUBDIRECTORIES)
    foreach(SUBDIR IN LISTS SUBDIRS)
        print_all_targets("${SUBDIR}")
    endforeach()
endfunction()

print_all_targets(.)

if (PYTHON_BINDINGS)
    add_subdirectory(python)
endif ()


add_executable(${EXE_NAME} main.cpp)
target_link_libraries(${EXE_NAME} ${LIB_NAME})
enable_lto(${EXE_NAME})
target_add_compilation_flags(${EXE_NAME})

file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data)

if (TESTS)
    message("COMPILING TESTS")
    target_compile_definitions(${LIB_NAME} PUBLIC TEST=42)
    add_subdirectory(tests EXCLUDE_FROM_ALL)
endif ()
