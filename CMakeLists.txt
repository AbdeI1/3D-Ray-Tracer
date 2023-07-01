cmake_minimum_required (VERSION 3.0)
project (3DRayTracer)

file(GLOB SOURCES
  src/*.cc
  main.cc
)

add_executable(raytracer ${SOURCES})

target_include_directories(raytracer
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)
