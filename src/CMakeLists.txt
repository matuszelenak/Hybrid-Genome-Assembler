cmake_minimum_required(VERSION 3.13)
project(src)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")

FIND_PACKAGE( Boost 1.65 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

add_executable(src kmer_categorizer.cpp SequenceReader.cpp SequenceReader.h KmerIterator.cpp KmerIterator.h consts.h structs.h KmerAnalysisWriter.cpp KmerAnalysisWriter.h ReadCategorizer.cpp ReadCategorizer.h)


TARGET_LINK_LIBRARIES( src LINK_PUBLIC ${Boost_LIBRARIES} stdc++fs pthread)