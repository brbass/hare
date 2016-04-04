cmake_minimum_required (VERSION 3.1)

project(Hare)

# compiler options
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_VERBOSE_MAKEFILE off)
set(CMAKE_BUILD_TYPE Debug)

# standard flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -o3")
set(CMAKE_INSTALL_PREFIX /home/brbass/research/hare/bin)

file(GLOB SOURCES "src/*.cc" "lib/*.cpp")

include_directories(src lib)

add_executable(hare ${SOURCES})

install(TARGETS hare DESTINATION .)

