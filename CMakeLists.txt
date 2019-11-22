# CMake version.
cmake_minimum_required(VERSION 3.10)

# Project.
project(layered-sqt VERSION 0.1 LANGUAGES CXX)

# Set release.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set compile flags.
set(CMAKE_CXX_FLAGS "-Wall -Wextra -march=native -mtune=native")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Require threads.
find_package(Threads REQUIRED)

# Glob sources.
file(GLOB_RECURSE SOURCES RELATIVE ${CMAKE_SOURCE_DIR} "src/*.cpp")

# Add executable.
add_executable(layered-sqt ${SOURCES})

# Set C++17.
set_target_properties(
    layered-sqt
    PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED True
    )

# Include directories.
target_include_directories(
    layered-sqt
    PUBLIC
    "${CMAKE_SOURCE_DIR}/include"
    "${CMAKE_SOURCE_DIR}/lib/preform/include"
    )

# Link threads.
target_link_libraries(
    layered-sqt ${CMAKE_THREAD_LIBS_INIT})
