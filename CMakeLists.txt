cmake_minimum_required(VERSION 2.6)
project(axi_psdf)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(CMAKE_CXX_FLAGS "-O3")
#set(CMAKE_CXX_FLAGS "-g")

file(GLOB srcs "src/*.cpp")
file(GLOB halos "src/halos/*.cpp")
file(GLOB baryons "src/baryons/*.cpp")

find_package(GSL REQUIRED)
add_executable(axi_psdf main.cpp ${halos} ${baryons} ${srcs})
target_link_libraries(axi_psdf GSL::gsl GSL::gslcblas pthread)

add_executable(tests tests.cpp ${halos} ${baryons} ${srcs})
target_link_libraries(tests GSL::gsl GSL::gslcblas pthread CppUTest)

install(TARGETS axi_psdf RUNTIME DESTINATION bin)
