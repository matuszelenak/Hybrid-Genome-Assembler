cmake_minimum_required(VERSION 3.13)
project(src)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -fconcepts")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -fconcepts")

FIND_PACKAGE( Boost 1.65 COMPONENTS program_options date_time system REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

add_subdirectory(lib/fmt EXCLUDE_FROM_ALL)
add_subdirectory(lib/robin-map EXCLUDE_FROM_ALL)


file(GLOB COMMON_SRC_FILES ${PROJECT_SOURCE_DIR}/common/*.cpp)
file(GLOB LIB_SRC_FILES ${PROJECT_SOURCE_DIR}/lib/*.*pp)

# Kmer occurrence executable
add_executable(occurrences kmer_occurrences.cpp ${COMMON_SRC_FILES} ${LIB_SRC_FILES})

TARGET_LINK_LIBRARIES( occurrences LINK_PUBLIC ${Boost_LIBRARIES} stdc++fs pthread fmt)
target_link_libraries( occurrences PRIVATE tsl::robin_map)



## Clustering executable
add_executable(clustering read_clustering.cpp ${COMMON_SRC_FILES} ${LIB_SRC_FILES})

TARGET_LINK_LIBRARIES( clustering LINK_PUBLIC ${Boost_LIBRARIES} stdc++fs pthread fmt)
target_link_libraries( clustering PRIVATE tsl::robin_map)


## Locating characteristic kmers in reference genome
add_executable(locate locate_in_reference.cpp ${COMMON_SRC_FILES} ${LIB_SRC_FILES})
TARGET_LINK_LIBRARIES( locate LINK_PUBLIC ${Boost_LIBRARIES} stdc++fs pthread fmt)
target_link_libraries( locate PRIVATE tsl::robin_map)


add_executable(subsample subsample.cpp common/SequenceRecordIterator.cpp common/Utils.cpp)
TARGET_LINK_LIBRARIES( subsample LINK_PUBLIC ${Boost_LIBRARIES} stdc++fs pthread fmt)