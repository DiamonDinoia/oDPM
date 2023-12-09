
find_package(OpenMP REQUIRED)
list(APPEND LIBRARIES OpenMP::OpenMP_CXX)
find_package(Threads REQUIRED)
list(APPEND LIBRARIES Threads::Threads)

if (OpenMP_CXX_FOUND)
    list(APPEND FLAGS ${OpenMP_CXX_FLAGS})
endif ()
