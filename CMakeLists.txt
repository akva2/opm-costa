cmake_minimum_required(VERSION 3.5)

project(opm-CoSTA)

cmake_policy(SET CMP0057 NEW)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(opm-common REQUIRED)
find_package(opm-grid REQUIRED)
find_package(opm-simulators REQUIRED)

set(CMAKE_CXX_STANDARD 17)

find_package(Python3 COMPONENTS Development Interpreter)
find_package(pybind11 REQUIRED)

pybind11_add_module(opm_CoSTA opm_CoSTA.cpp)
target_sources(opm_CoSTA PRIVATE flow_onephase_CoSTA.cpp
                                 moduleVersion.cpp 
                                 CoSTAModule.hpp)

target_include_directories(opm_CoSTA PRIVATE ${opm-simulators_INCLUDE_DIRS})

target_link_libraries(opm_CoSTA PUBLIC opmcommon opmgrid opmsimulators)
