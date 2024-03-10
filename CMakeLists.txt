cmake_minimum_required(VERSION 3.27)
project(NM_sopr_grad_laba)

set(CMAKE_CXX_STANDARD 17)

set(Boost_INCLUDE_DIR c:/MinGW_with_boost/MinGW/include/)
set(Boost_LIBRARY_DIR c:/MinGW_with_boost/MinGW/lib/)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

add_executable(NM_sopr_grad_laba main.cpp)

target_link_libraries(NM_sopr_grad_laba ${Boost_LIBRARY_DIR})
