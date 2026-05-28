cmake_minimum_required(VERSION 3.15)

project(BOSRayTracePreview LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(bos_ray_trace main.cpp)

if (MSVC)
    target_compile_options(bos_ray_trace PRIVATE /utf-8)
endif()
