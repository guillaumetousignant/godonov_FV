cmake_minimum_required(VERSION 3.12)
# 3.18 enables CUDA architecture flags

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "Do not build in-source. Please remove CMakeCache.txt and the CMakeFiles/ directory. Then build out-of-source.")
endif()

project(assignment_3_project VERSION 1.0.0 LANGUAGES C CXX)

add_executable(assignment_3 
    src/main.cpp)

target_include_directories(assignment_3 
    PRIVATE
        include/
)

target_compile_features(assignment_3 PRIVATE cxx_std_20)

option(BUILD_TESTING "Build all tests." OFF)

if (BUILD_TESTING)
    include(FetchContent)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v2.13.3
    )
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/contrib)

    include(CTest)
    include(Catch)

    add_executable(unit_tests 
        tests/test_main.cpp
        tests/example_test.cpp)
    target_link_libraries(unit_tests PRIVATE 
        Catch2::Catch2)
    target_compile_features(unit_tests PRIVATE cxx_std_20)

    catch_discover_tests(unit_tests)
endif()

install(TARGETS assignment_3 DESTINATION bin)