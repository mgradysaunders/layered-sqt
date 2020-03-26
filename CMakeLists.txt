# CMake version.
cmake_minimum_required(VERSION 3.10)

# Project.
project(layered-sqt VERSION 0.1 LANGUAGES C CXX)

# Set C11/C++17.
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set C++ filesystem library.
if(NOT CXX_FILESYSTEM)
    set(CXX_FILESYSTEM stdc++fs)
    message(STATUS "Setting CXX_FILESYSTEM=stdc++fs")
endif()

# Set release.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    message(STATUS "Setting CMAKE_BUILD_TYPE=Release")
endif()

# Set C++ compile flags.
set(CMAKE_CXX_FLAGS "-Wall -Wextra -march=native -mtune=native")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")

# Set C compile flags.
set(CMAKE_C_FLAGS "-Wall -Wextra -march=native -mtune=native")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG")

# Ensure CXX compiler supports CXX filesystem.
try_compile(
    HAS_CXX_FILESYSTEM
    "${CMAKE_BINARY_DIR}/tmp"
    "${CMAKE_SOURCE_DIR}/tests/has_filesystem.cpp"
    CMAKE_FLAGS 
    -DCMAKE_CXX_STANDARD=17 
    -DCMAKE_CXX_STANDARD_REQUIRED=ON
    LINK_LIBRARIES ${CXX_FILESYSTEM}
    )
if(HAS_CXX_FILESYSTEM)
    message(
        STATUS 
        "${CMAKE_CXX_COMPILER} supports C++17 <filesystem>")
else()
    message(
        FATAL_ERROR
        "${CMAKE_CXX_COMPILER} does not support C++17 <filesystem>")
endif()

# Require threads.
find_package(Threads REQUIRED)

# Add library.
add_library(
    layered-sqt-core 
    STATIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/file_data.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/layer.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/layer/lambertian.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/layer/microsurface_conductive.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/layer/microsurface_dielectric.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/layer/microsurface_lambertian.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/layer/null.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/layer/oren_nayar_diffuse.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/layered_assembly.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/medium.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/medium/henyey_greenstein.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/medium/henyey_greenstein2.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/medium/rayleigh.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/medium/sggx.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/progress_bar.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/rrss.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/tri.cpp"
    )
target_include_directories(
    layered-sqt-core
    PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/preform/include"
    )

# Add executable.
add_executable(
    layered-sqt
    "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp"
    )
target_include_directories(
    layered-sqt
    PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/preform/include"
    )
target_link_libraries(layered-sqt layered-sqt-core)
target_link_libraries(layered-sqt ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(layered-sqt ${CXX_FILESYSTEM})

# Add executable.
add_executable(
    layered-sqt-lssinfo
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lssinfo.cpp"
    )
target_include_directories(
    layered-sqt-lssinfo
    PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/preform/include"
    )
target_link_libraries(layered-sqt-lssinfo layered-sqt-core)

# Add executable.
add_executable(
   layered-sqt-lssview
    "${CMAKE_CURRENT_SOURCE_DIR}/src/stb_image_write.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lssview.cpp"
    )
target_include_directories(
    layered-sqt-lssview
    PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/preform/include"
    )
target_link_libraries(layered-sqt-lssview layered-sqt-core)
