set(FLAGS_DEBUG -O0 -fstack-protector-all -Wall -Wextra -pedantic -fno-inline -fno-omit-frame-pointer -g)
set(FLAGS_RELEASE -O3 -ffp-contract=fast -freciprocal-math
        -fmerge-all-constants -fno-math-errno -march=native -funroll-loops
        -ftree-vectorize -fno-trapping-math -fassociative-math -ffinite-math-only
        -fno-signed-zeros)
set(FLAGS_RelWithDebInfo ${FLAGS_RELEASE} -g -fno-omit-frame-pointer)

if (NOT (CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
    list(APPEND FLAGS_RELEASE -fcx-limited-range)
endif ()

list(APPEND FLAGS $<$<CONFIG:DEBUG>:${FLAGS_DEBUG}>
        $<$<CONFIG:Release>:${FLAGS_RELEASE}>
        $<$<CONFIG:RelWithDebInfo>:${FLAGS_RelWithDebInfo}>)

function(target_add_compilation_flags target)
    set(old_flags ${FLAGS})
    message(STATUS "COMPILATION FLAGS ${old_flags}")
    target_compile_options(${target} PUBLIC  "$<$<COMPILE_LANGUAGE:CXX>:${old_flags}>")
    if (GPU)
        if (NOT "${old_flags}" STREQUAL "")
            string(REPLACE ";" "," CUDA_flags "${old_flags}")
            string(REPLACE "-pedantic" "-Wno-pedantic" CUDA_flags "${CUDA_flags}")
            target_compile_options(${target} PUBLIC "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${CUDA_flags}>")
        endif ()
    endif ()
endfunction()