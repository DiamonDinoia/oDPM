file(GLOB_RECURSE VALIDATION_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

add_library(DPMValidation STATIC ${VALIDATION_SOURCES})
target_include_directories(DPMValidation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(DPMValidation PUBLIC ${LIB_NAME})
