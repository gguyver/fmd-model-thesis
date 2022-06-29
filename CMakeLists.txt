cmake_minimum_required(VERSION 3.14)
project(model_fmd)

include_directories(include external/yaml-cpp/include/yaml-cpp)

add_executable(model_fmd main.cpp)
set_target_properties(model_fmd PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
target_compile_options(model_fmd PRIVATE -Werror -Wall -Wextra -pedantic-errors -Wconversion -Wsign-conversion)

add_subdirectory(src)
target_link_libraries(model_fmd PRIVATE sources)


set_property(GLOBAL PROPERTY CTEST_TARGETS_ADDED 1) # Remove CTest Configurations from CLion
add_subdirectory(external/yaml-cpp)
target_link_libraries(model_fmd PRIVATE yaml-cpp)
target_link_libraries(sources PRIVATE yaml-cpp)

iF (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    enable_testing()
    include(GoogleTest)
    add_subdirectory(tests)
endif ()

