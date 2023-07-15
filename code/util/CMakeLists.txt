cmake_minimum_required(VERSION 3.23)
project(bsax)

set(CMAKE_CXX_STANDARD 11)

file(GLOB UTIL_LIB_SRC   # 
        ./*.cpp
        ./*.h
        )


# 
add_library(util SHARED ${UTIL_LIB_SRC})