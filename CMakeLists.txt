cmake_minimum_required(VERSION 2.6)
project(axi_psdf)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_SHARED_MODULE_PREFIX "")
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS "-O3")
#set(CMAKE_CXX_FLAGS "-g")
set(PYBIND11_PYTHON_VERSION 2.7)

file(GLOB srcs "src/*.cpp")
file(GLOB halos "src/halos/*.cpp")
file(GLOB baryons "src/baryons/*.cpp")

find_package(GSL REQUIRED)

add_executable(aim main.cpp ${halos} ${baryons} ${srcs})
target_link_libraries(aim GSL::gsl GSL::gslcblas pthread)
install(TARGETS aim RUNTIME DESTINATION bin)

add_executable(tests tests.cpp ${halos} ${baryons} ${srcs})
target_link_libraries(tests GSL::gsl GSL::gslcblas pthread CppUTest)

find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(AIM wrapper.cpp ${halos} ${baryons} ${srcs})
target_link_libraries(AIM GSL::gsl GSL::gslcblas pthread)

