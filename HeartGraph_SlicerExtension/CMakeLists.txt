cmake_minimum_required(VERSION 2.8.7)
set(MODULE_NAME GraphCutSeg)

# Use Slicer
FIND_PACKAGE(Slicer REQUIRED)
INCLUDE(${Slicer_USE_FILE})

SET(CMAKE_BUILD_TYPE Debug)
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

# SlicerExecutionModel
FIND_PACKAGE(GenerateCLP REQUIRED)
INCLUDE(${GenerateCLP_USE_FILE})

SET(mesh2graph2surfVnet_SOURCE mesh2graph2surfVnet.cxx)
GENERATECLP(mesh2graph2surfVnet_SOURCE mesh2graph2surfVnet.xml)

ADD_EXECUTABLE(mesh2graph2surfVnet mesh2graph2surfVnet.cxx ResampleVolume.cxx graph.cpp maxflow.cpp)
TARGET_LINK_LIBRARIES(mesh2graph2surfVnet ${ITK_LIBRARIES})

ADD_LIBRARY(mesh2graph2surfVnetLib SHARED mesh2graph2surfVnet mesh2graph2surfVnet.cxx ResampleVolume.cxx graph.cpp maxflow.cpp)
SET_TARGET_PROPERTIES(mesh2graph2surfVnetLib PROPERTIES COMPILE_FLAGS "-Dmain=ModuleEntryPoint" )
TARGET_LINK_LIBRARIES(mesh2graph2surfVnetLib ${ITK_LIBRARIES})
