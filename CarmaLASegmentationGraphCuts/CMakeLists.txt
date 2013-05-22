cmake_minimum_required(VERSION 2.8.7)

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

set(MODULE_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  )

#-----------------------------------------------------------------------------
# LA Segmentation Using Graph Cuts Module
#-----------------------------------------------------------------------------
set(MODULE_NAME CarmaLASegmentationGraphCuts)
set(MODULE_ADDITIONAL_SRCS ResampleVolume.cxx graph.cpp maxflow.cpp)
#-----------------------------------------------------------------------------

SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  ADDITIONAL_SRCS ${MODULE_ADDITIONAL_SRCS}
  LOGO_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/../NamicCarmaLogo.h
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  )
