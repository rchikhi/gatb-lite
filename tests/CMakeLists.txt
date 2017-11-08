cmake_minimum_required(VERSION 2.8)

#project(gatbl-tests)
include_directories( "../thirdparty/catch/single_include" )
add_executable(tests "main.cpp" "kmer.cpp")
target_link_libraries(tests  -lrt -lpthread)
