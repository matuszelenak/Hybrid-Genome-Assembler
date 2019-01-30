cmake_minimum_required(VERSION 3.13)
project(src)

set(CMAKE_CXX_STANDARD 14)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")

FIND_PACKAGE( Boost 1.65 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

add_executable(src kmer_categorizer.cpp)


TARGET_LINK_LIBRARIES( src LINK_PUBLIC ${Boost_LIBRARIES} )