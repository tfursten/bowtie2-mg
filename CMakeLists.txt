cmake_minimum_required(VERSION 2.8)
project(bowtie2_mg)

add_definitions(-DWRAPPER_VERSION="1.0")
add_definitions(-DBUFFER_SIZE=1024)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(BEFORE include)

add_executable(bowtie2_mg_aligner aligner.cpp)
add_executable(bowtie2_mg_indexer indexer.cpp)
