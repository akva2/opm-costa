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
target_sources(opm_CoSTA PRIVATE flow_blackoil_CoSTA.cpp
                                 flow_gaswater_CoSTA.cpp
                                 flow_onephase_CoSTA.cpp
                                 ${opm-simulators_SOURCE_DIR}/opm/simulators/utils/moduleVersion.cpp
                                 CoSTAModule.hpp)

target_include_directories(opm_CoSTA PRIVATE ${opm-simulators_INCLUDE_DIRS} ${opm-simulators_DIR})

target_link_libraries(opm_CoSTA PUBLIC opmsimulators opmgrid opmcommon)

# OPM libraries are built with OpenMP when available. HAVE_OPENMP changes
# class layouts, so it must match or we get ODR violations and crashes.
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(opm_CoSTA PUBLIC OpenMP::OpenMP_CXX)
  target_compile_definitions(opm_CoSTA PRIVATE HAVE_OPENMP=1)
endif()

# OPM builds its libraries with asserts enabled (-UNDEBUG) also in release
# mode, and some classes change layout with NDEBUG. Match it.
target_compile_options(opm_CoSTA PRIVATE $<$<CONFIG:Release,RelWithDebInfo>:-UNDEBUG>)
