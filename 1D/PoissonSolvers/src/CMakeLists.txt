cmake_minimum_required(VERSION 3.22)

project(Poisson-1D VERSION 1.0)

add_executable(fd fd/sparse.cpp)

set(CMAKE_CXX_STANDARD 17)
