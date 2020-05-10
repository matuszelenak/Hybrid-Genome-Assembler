cmake_minimum_required(VERSION 3.13)
project(src)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -fconcepts")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -fconcepts")

FIND_PACKAGE( Boost 1.65 COMPONENTS program_options date_time system REQUIRED regex )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

add_subdirectory(lib/fmt EXCLUDE_FROM_ALL)

file(GLOB COMMON_SRC_FILES ${PROJECT_SOURCE_DIR}/common/*.cpp)
file(GLOB LIB_SRC_FILES ${PROJECT_SOURCE_DIR}/lib/*.*pp)
file(GLOB OCCURRENCES_SRC_FILES ${PROJECT_SOURCE_DIR}/occurrences/*.*pp)
file(GLOB CLUSTERING_SRC_FILES ${PROJECT_SOURCE_DIR}/clustering/*.*pp)

# Kmer occurrence executable
add_executable(occurrences kmer_occurrences.cpp ${OCCURRENCES_SRC_FILES} ${COMMON_SRC_FILES} ${LIB_SRC_FILES})
TARGET_LINK_LIBRARIES( occurrences LINK_PUBLIC ${Boost_LIBRARIES} stdc++fs pthread fmt)


## Clustering executable
add_executable(clustering read_clustering.cpp ${CLUSTERING_SRC_FILES} ${COMMON_SRC_FILES} ${LIB_SRC_FILES})
TARGET_LINK_LIBRARIES( clustering LINK_PUBLIC ${Boost_LIBRARIES} stdc++fs pthread fmt)


## Subsampling reads
add_executable(subsample subsample.cpp common/SequenceRecordIterator.cpp common/Utils.cpp)
TARGET_LINK_LIBRARIES( subsample LINK_PUBLIC ${Boost_LIBRARIES} stdc++fs pthread fmt)


## Locating discriminative kmers in reference
add_executable(find_in_reference find_in_reference.cpp ${COMMON_SRC_FILES} ${LIB_SRC_FILES})
TARGET_LINK_LIBRARIES( find_in_reference LINK_PUBLIC ${Boost_LIBRARIES} stdc++fs pthread fmt)


add_executable(jf_occurrences jellyfish_occurrences.cpp ${OCCURRENCES_SRC_FILES} ${COMMON_SRC_FILES} ${LIB_SRC_FILES})
TARGET_LINK_LIBRARIES( jf_occurrences LINK_PUBLIC ${Boost_LIBRARIES} stdc++fs pthread fmt)