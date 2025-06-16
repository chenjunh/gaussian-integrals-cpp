cmake_minimum_required(VERSION 3.10)
project(gaussian_integrals_cpp)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)

add_executable(test_hermite tests/test_hermite.cpp src/hermite.cpp)

enable_testing()
add_test(NAME HermiteTest COMMAND test_hermite)
