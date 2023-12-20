#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(dicomReader)

#----------------------------------------------------------------------------

set(DICOM_READER_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/include)
file(GLOB DICOM_READER_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
add_library(g4DicomReader STATIC ${DICOM_READER_SOURCES})
target_include_directories(g4DicomReader PUBLIC ${DICOM_READER_INCLUDE})
target_include_directories(g4DicomReader PRIVATE
        ${dcmtk_SOURCE_DIR}/dcmdata/include
        ${dcmtk_SOURCE_DIR}/ofstd/include
        ${dcmtk_SOURCE_DIR}/oflog/include
        ${DCMTK_BINARY_DIR}/config/include/
)

target_link_libraries(g4DicomReader PUBLIC DCMTK config G4run-static G4zlib-static dcmdata
        G4geomtext-static)
