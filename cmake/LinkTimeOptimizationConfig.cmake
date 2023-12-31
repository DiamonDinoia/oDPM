
include(CheckIPOSupported)
check_ipo_supported(RESULT supported OUTPUT error)
if (supported)
    if (POLICY CMP0069)
        cmake_policy(SET CMP0069 NEW)
        set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
    endif ()
    message(STATUS "IPO / LTO enabled")
else ()
    message(STATUS "IPO / LTO not supported: <${error}>")
endif ()

function(enable_lto target)
    if (supported)
        set_target_properties(${target} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${ENABLE_LTO})
    else ()
        message(WARNING "IPO / LTO not supported: <${error}>")
    endif ()
endfunction()