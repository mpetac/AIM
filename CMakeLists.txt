cmake_minimum_required(VERSION 2.6)
project(axi_psdf)

add_executable(axi_psdf main.cpp)

install(TARGETS axi_psdf RUNTIME DESTINATION bin)
