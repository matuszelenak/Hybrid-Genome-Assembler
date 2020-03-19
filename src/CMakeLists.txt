cmake_minimum_required(VERSION 3.13)
project(src)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")

FIND_PACKAGE( Boost 1.65 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

add_subdirectory(fmt EXCLUDE_FROM_ALL)

# add_executable(src kmer_categorizer.cpp SequenceReader.cpp SequenceReader.h KmerIterator.cpp KmerIterator.h consts.h structs.h KmerAnalysisWriter.cpp KmerAnalysisWriter.h ReadCategorizer.cpp ReadCategorizer.h)

add_executable(occurrences KmerOccurrences.cpp SequenceReader.cpp SequenceReader.h KmerIterator.cpp KmerIterator.h Utils.cpp Utils.h structs.h consts.h)

add_executable(categorizer CategorizeReads.cpp GenomeReadCluster.h GenomeReadCluster.cpp)

TARGET_LINK_LIBRARIES( occurrences LINK_PUBLIC ${Boost_LIBRARIES} stdc++fs pthread jsoncpp fmt)
