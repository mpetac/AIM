cmake_minimum_required(VERSION 2.6)
project(axi_psdf)

add_executable(axi_psdf main.cpp Baryons_A.cpp Parametric_funcs.cpp)
add_executable(tests Tests.cpp Parametric_funcs.cpp)
target_link_libraries(tests CppUTest )

install(TARGETS axi_psdf RUNTIME DESTINATION bin)
