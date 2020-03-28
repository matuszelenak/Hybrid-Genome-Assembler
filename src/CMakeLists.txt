cmake_minimum_required(VERSION 3.13)
project(src)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

FIND_PACKAGE( Boost 1.65 COMPONENTS program_options date_time system REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

add_subdirectory(fmt EXCLUDE_FROM_ALL)
add_subdirectory(robin-map EXCLUDE_FROM_ALL)

add_executable(pipeline MainPipeline.cpp KmerAnalysis.cpp KmerIterator.cpp Utils.cpp DNAStructures.cpp SequenceRecordIterator.cpp ReadClusteringEngine.cpp)

TARGET_LINK_LIBRARIES( pipeline LINK_PUBLIC ${Boost_LIBRARIES} stdc++fs pthread fmt)
target_link_libraries( pipeline PRIVATE tsl::robin_map)

add_executable(playground Playground.cpp KmerCountingBloomFilter.cpp smhasher/MurmurHash3.cpp KmerAnalysis.cpp KmerIterator.cpp SequenceRecordIterator.cpp Utils.cpp)
target_link_libraries( playground LINK_PUBLIC ${Boost_LIBRARIES} fmt pthread)
target_link_libraries( playground PRIVATE tsl::robin_map)
