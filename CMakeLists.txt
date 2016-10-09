cmake_minimum_required(VERSION 3.6)
project(lab_4_5)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -fopenmp")

set(SOURCE_FILES main.cpp Sequential.cpp test.cpp OptimizedParallelMultiplication.cpp DefaultParallel.cpp Assignment_120339P_120168N.cpp)
add_executable(lab_4_5 ${SOURCE_FILES})