cmake_minimum_required(VERSION 3.23)
project(bsax)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -mavx2")

file(GLOB BTREE_LIB_SRC   # 
        ./*.cpp
        ./*.h)

include_directories(../sax/include)
include_directories(../util)

# 
add_library(btree SHARED ${BTREE_LIB_SRC})  # esbtree.so,SHARED,