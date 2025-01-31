cmake_minimum_required(VERSION 3.4)
project(BBS CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")


include(FetchContent)

# seqan3 for the basic skeleton for the mapper
FetchContent_Declare(seqan3
    GIT_REPOSITORY https://github.com/seqan/seqan3.git
    GIT_TAG        3.3.0
)
FetchContent_MakeAvailable(seqan3)

# Sharg for argument parsing
FetchContent_Declare(sharg
    GIT_REPOSITORY https://github.com/seqan/sharg-parser.git
    GIT_TAG        1.1.1
)
FetchContent_MakeAvailable(sharg)

add_executable(bbs src/main.cpp src/bbs.hpp src/utils.hpp)
target_link_libraries(bbs PUBLIC seqan3::seqan3 sharg::sharg)