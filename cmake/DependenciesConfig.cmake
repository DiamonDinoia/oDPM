
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    set(CMAKE_CUDA_COMPILER_LAUNCHER "${CCACHE_PROGRAM}") # CMake 3.9+
endif()

set(MOLD_PROGRAM mold)
if (WIN32)
    set(MOLD_PROGRAM mold.exe)
endif ()

find_program(MOLD_PROGRAM ${MOLD_PROGRAM})
if(MOLD_PROGRAM)
    set(CMAKE_LINKER_TYPE MOLD)
endif ()

CPMAddPackage(
        NAME mdspan
        GIT_REPOSITORY https://github.com/kokkos/mdspan.git
        GIT_TAG mdspan-0.6.0
        VERSION 0.6.0
        GIT_SHALLOW YES
        GIT_PROGRESS YES
        EXCLUDE_FROM_ALL YES
        SYSTEM
        OPTIONS
        "MDSPAN_ENABLE_CUDA ${GPU}"

)
target_compile_definitions(mdspan INTERFACE "MDSPAN_USE_PAREN_OPERATOR=1")
list(APPEND LIBRARIES mdspan)

CPMAddPackage(
        NAME abseil
        GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
        GIT_TAG 20230802.1
        GIT_PROGRESS YES
        GIT_SHALLOW YES
        EXCLUDE_FROM_ALL YES
        SYSTEM
        OPTIONS "ABSL_PROPAGATE_CXX_STD ON"
                "ABSL_BUILD_TESTING OFF"
                "ABSL_BUILD_TEST_HELPERS OFF"
)

CPMAddPackage(
        NAME dcmtk
        GIT_REPOSITORY https://github.com/DCMTK/dcmtk.git
        GIT_TAG DCMTK-3.6.7
        GIT_PROGRESS YES
        GIT_SHALLOW YES
        EXCLUDE_FROM_ALL YES
        SYSTEM
        OPTIONS
        "DCMTK_DEFAULT_DICT builtin"
)

if (NOT GEANT4_INSTALL_DATADIR)
    message(STATUS "GEANT4_INSTALL_DATADIR not set, using default: ${PROJECT_SOURCE_DIR}/geant4")
    set(GEANT4_INSTALL_DATADIR ${PROJECT_SOURCE_DIR}/geant4)
endif ()

if (NOT EXISTS ${GEANT4_INSTALL_DATADIR})
    file(MAKE_DIRECTORY ${GEANT4_INSTALL_DATADIR})
endif ()

list(APPEND DATASETS G4NDL.4.7.tar.gz G4EMLOW.8.2.tar.gz G4PhotonEvaporation.5.7.tar.gz G4RadioactiveDecay.5.6.tar.gz
        G4PARTICLEXS.4.0.tar.gz G4RealSurface.2.2.tar.gz G4SAIDDATA.2.0.tar.gz G4ABLA.3.1.tar.gz G4INCL.1.0.tar.gz
        G4ENSDFSTATE.2.3.tar.gz G4TENDL.1.4.tar.gz G4PII.1.3.tar.gz)

foreach (dataset ${DATASETS})
    # strip version and extension from dataset name
    string(REPLACE "." ";" dataset_name ${dataset})
    list(GET dataset_name 0 dataset_name)
    if (NOT EXISTS ${GEANT4_INSTALL_DATADIR}/${dataset_name})
        file(MAKE_DIRECTORY ${GEANT4_INSTALL_DATADIR}/${dataset_name})
        message(STATUS "Downloading ${dataset_name}... ${dataset} into ${GEANT4_INSTALL_DATADIR}/${dataset_name}")
        CPMAddPackage(
                NAME ${dataset_name}
                URL https://cern.ch/geant4-data/datasets/${dataset}
                EXCLUDE_FROM_ALL YES YES
                DOWNLOAD_ONLY YES
        )
        file(COPY ${${dataset_name}_SOURCE_DIR} DESTINATION ${GEANT4_INSTALL_DATADIR}/)
        file(RENAME ${${dataset_name}_SOURCE_DIR} ${GEANT4_INSTALL_DATADIR}/${dataset_name})
    endif ()
endforeach ()

CPMAddPackage(
        NAME geant4
        GIT_REPOSITORY https://github.com/Geant4/geant4.git
        GIT_TAG v11.2.0
        GIT_SHALLOW YES
        GIT_PROGRESS YES
        EXCLUDE_FROM_ALL YES
        SYSTEM
        OPTIONS
        "GEANT4_USE_SYSTEM_EXPAT OFF"
        "GEANT4_USE_SYSTEM_ZLIB OFF"
        "BUILD_SHARED_LIBS OFF"
        "BUILD_STATIC_LIBS ON"
        "GEANT4_INSTALL_DATA OFF"
        "GEANT4_INSTALL_DATASETS_TENDL OFF"
        "GEANT4_INSTALL_DATADIR ${GEANT4_INSTALL_DATADIR}"
        "GEANT4_USE_QT OFF"
)

CPMAddPackage(
        NAME dpmg4cpp
        GIT_REPOSITORY git@github.com:DiamonDinoia/dpm-g4cpp.git
        GIT_TAG master
        GIT_PROGRESS YES
        GIT_SHALLOW YES
        EXCLUDE_FROM_ALL YES
        SYSTEM
)

CPMAddPackage(
        NAME mixmax
        GIT_REPOSITORY https://github.com/DiamonDinoia/mixmaxCUDA.git
        GIT_TAG main
        GIT_PROGRESS YES
        GIT_SHALLOW YES
        EXCLUDE_FROM_ALL YES
        SYSTEM
)

CPMAddPackage(
        NAME Boost
        GIT_REPOSITORY https://github.com/boostorg/boost.git
        GIT_TAG boost-1.83.0
        GIT_PROGRESS YES
        GIT_SHALLOW YES
        EXCLUDE_FROM_ALL YES
        SYSTEM
)

CPMAddPackage(
        NAME yagit
        GIT_REPOSITORY https://github.com/DiamonDinoia/yagit.git
        GIT_TAG no-dependency
        GIT_PROGRESS YES
        GIT_SHALLOW YES
        EXCLUDE_FROM_ALL YES
        SYSTEM
        DOWNLOAD_ONLY YES
        OPTIONS
        "BUILD_EXAMPLES NO"
        "BUILD_TESTING NO"
        "BUILD_PERFORMANCE_TESTING NO"
)

add_subdirectory(${yagit_SOURCE_DIR}/gi_core ${yagit_BINARY_DIR} EXCLUDE_FROM_ALL)

list(APPEND LIBRARIES absl::inlined_vector absl::span absl::optional mixmax dpm_data g4Dicom g4DicomReader
        G4physicslists-static G4vis_management-static G4Tree-static G4FR-static G4visHepRep-static G4RayTracer-static
        G4VRML-static G4GMocren-static G4ToolsSG-static G4interfaces-static gi_core
)
