cmake_minimum_required(VERSION 4.0)
project(group7)

set(CMAKE_CXX_STANDARD 20)

add_executable(group7 src/main.cpp
        include/graph/Graph.h
        include/graph/Edge.h
        src/graph/Graph.cpp
        include/util/Timer.h
)

target_include_directories(group7 PRIVATE include)
include_directories(mtx)

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
        GraphTest
        tests/GraphTest.cpp
        src/graph/Graph.cpp
)
target_link_libraries(
        GraphTest
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(GraphTest)
