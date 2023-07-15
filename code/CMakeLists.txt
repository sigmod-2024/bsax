cmake_minimum_required(VERSION 3.23)
project(bsax)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -mavx2")

include_directories(./sax/include)
include_directories(./util)

add_subdirectory(sax)
add_subdirectory(util)
add_subdirectory(btree)
add_subdirectory(i_binary)
add_subdirectory(b_binary)
add_subdirectory(sigtree)


add_executable(walk_through test_walk_through.cpp)
target_link_libraries(walk_through sax util)

add_executable(test_b_btree bak/test_btree.cpp)
target_link_libraries(test_b_btree btree sax util)

add_executable(test_binary bak/test_binary.cpp)
target_link_libraries(test_binary b_binary sax util)


#add_subdirectory(bak/i_btree)
#add_executable(test_i_btree bak/test_i_btree.cpp)
#target_link_libraries(test_i_btree bak/i_btree sax util)


add_executable(exp2_btree exp2_btree.cpp)
target_link_libraries(exp2_btree btree sax util)
add_executable(exp2_binary exp2_binary.cpp)
target_link_libraries(exp2_binary i_binary b_binary sax util)
add_executable(exp3 exp3.cpp)
target_link_libraries(exp3 i_binary b_binary btree sax util)
add_executable(exp4 exp4.cpp)
target_link_libraries(exp4 i_binary b_binary btree sax util)
#add_executable(exp5_btree bak/exp5_btree.cpp)
#target_link_libraries(exp5_btree btree sax util)
#add_executable(exp5_binary bak/exp5_binary.cpp)
#target_link_libraries(exp5_binary i_binary b_binary sax util)
#add_executable(exp5_sigtree bak/exp5_sigtree.cpp)
#target_link_libraries(exp5_sigtree sigtree sax util)
add_executable(exp5 exp5.cpp)
target_link_libraries(exp5 i_binary b_binary btree sigtree sax util)

add_executable(generate_sax generate_sax.cpp)
target_link_libraries(generate_sax sax util)
add_executable(exp2_histogram exp2_histogram.cpp)
add_executable(exp4_test exp4_test.cpp)
target_link_libraries(exp4_test sax util)

add_executable(reorder reorder.cpp)
target_link_libraries(reorder sax)

add_executable(test_sax test_sax.cpp)
target_link_libraries(test_sax sax)

add_executable(generate_ts generate_ts.cpp)
target_link_libraries(generate_ts sax util)

add_executable(find_read find_read.cpp)
target_link_libraries(find_read sax util)

add_executable(exp4_test1 exp4_test1.cpp)
target_link_libraries(exp4_test1 sax util)

#add_executable(main main.cpp)
#target_link_libraries(main sax)