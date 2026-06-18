cmake_minimum_required (VERSION 3.20...3.31)
project (BABjarneJehn CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable (get_distribution.cpp)


add_subdirectory(rshash)

