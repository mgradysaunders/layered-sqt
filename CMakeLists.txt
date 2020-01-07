# CMake version.
cmake_minimum_required(VERSION 3.10)

# Project.
project(layered-sqt VERSION 0.1 LANGUAGES C CXX)

# Set release.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set C compile flags.
set(CMAKE_C_FLAGS "-Wall -Wextra -march=native -mtune=native")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG")

# Set C++ compile flags.
set(CMAKE_CXX_FLAGS "-Wall -Wextra -march=native -mtune=native")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Require threads.
find_package(Threads REQUIRED)

# Set C++17.
function(set_target_cxx17 TARGET_ARG)
    set_target_properties(
        ${TARGET_ARG}
        PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED True
        )
endfunction(set_target_cxx17)

# Set common include directories.
function(set_target_common_include_directories TARGET_ARG)
    target_include_directories(
        ${TARGET_ARG}
        PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/lib/preform/include"
        )
endfunction(set_target_common_include_directories)

# Core sources.
set(CORE_SOURCES
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
    "${CMAKE_CURRENT_SOURCE_DIR}/src/medium/rayleigh.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/medium/sggx.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/progress_bar.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/rrss.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/tri.cpp"
    )

# Add library.
add_library(stbiw STATIC "${CMAKE_CURRENT_SOURCE_DIR}/src/stb_image_write.c")
set_target_properties(
    stbiw
    PROPERTIES
    C_STANDARD 11
    C_STANDARD_REQUIRED True
    )

# Add library.
add_library(layered-sqt-core STATIC ${CORE_SOURCES})
set_target_cxx17(layered-sqt-core)
set_target_common_include_directories(layered-sqt-core)

# Add executable.
add_executable(
    layered-sqt
    "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp"
    )
set_target_cxx17(layered-sqt)
set_target_common_include_directories(layered-sqt)
target_link_libraries(layered-sqt layered-sqt-core)
target_link_libraries(layered-sqt ${CMAKE_THREAD_LIBS_INIT})

# Add executable.
add_executable(
    layered-sqt-lssinfo
    "${CMAKE_CURRENT_SOURCE_DIR}/src/main_lssinfo.cpp"
    )
set_target_cxx17(layered-sqt-lssinfo)
set_target_common_include_directories(layered-sqt-lssinfo)
target_link_libraries(layered-sqt-lssinfo layered-sqt-core)

# Add executable.
add_executable(
    layered-sqt-lssview
    "${CMAKE_CURRENT_SOURCE_DIR}/src/main_lssview.cpp"
    )
set_target_cxx17(layered-sqt-lssview)
set_target_common_include_directories(layered-sqt-lssview)
target_link_libraries(layered-sqt-lssview layered-sqt-core)
target_link_libraries(layered-sqt-lssview stbiw)
