include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(basic_tests main_basic_tests.cpp node_tests.cpp config_tests.cpp common_tests.cpp matrix_tests.cpp report_tests.cpp cell_tests.cpp grid_tests.cpp point_tests.cpp disease_compartment_tests.cpp serotype_compartmental_model_tests.cpp)
set_target_properties(basic_tests PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)

target_link_libraries(basic_tests gtest gtest_main)
target_link_libraries(basic_tests sources)
target_link_libraries(basic_tests yaml-cpp)
