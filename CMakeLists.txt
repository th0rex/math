cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(math LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(math INTERFACE)
target_include_directories(math INTERFACE include/)

add_executable(dummy src/dummy.cpp)
target_link_libraries(dummy math)

add_executable(16_2_elgamal src/16_2_elgamal.cpp)
target_link_libraries(16_2_elgamal math)

add_executable(25_1 src/25_1.cpp)
target_link_libraries(25_1 math)
