# using the highest C++ standard supported by the compiler
# at least c++20 is required

set(ODPM_CXX_STANDARD DETECT CACHE STRING "Override the default CXX_STANDARD to compile with.")
set_property(CACHE ODPM_CXX_STANDARD PROPERTY STRINGS DETECT)

if(ODPM_CXX_STANDARD STREQUAL "DETECT")
    foreach(CXX_STANDARD_VAR 23;20)
        if("cxx_std_${CXX_STANDARD_VAR}" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
            message(STATUS "Detected support for C++${CXX_STANDARD_VAR} standard")
            set(ODPM_CXX_STANDARD ${CXX_STANDARD_VAR})
            break()
        endif()
    endforeach()
    if (NOT ODPM_CXX_STANDARD)
        message(FATAL_ERROR "No supported C++ standard available")
    endif ()
    message(STATUS "ODPM Using C++${ODPM_CXX_STANDARD} standard")
endif()
