cmake_minimum_required(VERSION 3.13)
project(src)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -fconcepts")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

FIND_PACKAGE( Boost 1.65 COMPONENTS program_options date_time system REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

add_subdirectory(fmt EXCLUDE_FROM_ALL)
add_subdirectory(robin-map EXCLUDE_FROM_ALL)


file(GLOB COMMON_SRC_FILES ${PROJECT_SOURCE_DIR}/common/*.cpp)
file(GLOB LIB_SRC_FILES ${PROJECT_SOURCE_DIR}/lib/*.*pp)

# Testing executable
file(GLOB TESTING_SRC_FILES ${PROJECT_SOURCE_DIR}/testing/*.cpp)
add_executable(pipeline_test ${COMMON_SRC_FILES} ${LIB_SRC_FILES} ${TESTING_SRC_FILES})

TARGET_LINK_LIBRARIES( pipeline_test LINK_PUBLIC ${Boost_LIBRARIES} stdc++fs pthread fmt)
target_link_libraries( pipeline_test PRIVATE tsl::robin_map)



## Production executable
file(GLOB PRODUCTION_SRC_FILES ${PROJECT_SOURCE_DIR}/production/*.cpp)
add_executable(pipeline ${COMMON_SRC_FILES} ${LIB_SRC_FILES} ${PRODUCTION_SRC_FILES})

TARGET_LINK_LIBRARIES( pipeline LINK_PUBLIC ${Boost_LIBRARIES} stdc++fs pthread fmt)
target_link_libraries( pipeline PRIVATE tsl::robin_map)




add_executable(playground playground/Playground.cpp ${COMMON_SRC_FILES} ${LIB_SRC_FILES})
TARGET_LINK_LIBRARIES( playground LINK_PUBLIC ${Boost_LIBRARIES} stdc++fs pthread fmt)
target_link_libraries( playground PRIVATE tsl::robin_map)