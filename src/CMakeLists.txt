include_directories(include)

add_library(sources Cell.cpp Common.cpp Config.cpp DiseaseCompartment.cpp Grid.cpp Node.cpp Point.cpp Reporting.cpp SerotypeCompartmentalModel.cpp)

set_target_properties(sources PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
target_compile_options(sources PRIVATE -Werror -Wall -Wextra -pedantic-errors -Wconversion -Wsign-conversion)