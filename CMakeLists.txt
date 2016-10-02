cmake_minimum_required(VERSION 3.6)
project(lab_4_5)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

set(SOURCE_FILES main.cpp parallelize_for_loop.cpp)
add_executable(lab_4_5 ${SOURCE_FILES})