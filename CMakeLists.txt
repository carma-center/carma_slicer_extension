cmake_minimum_required(VERSION 2.8.7)

#-----------------------------------------------------------------------------
set(EXTENSION_NAME carma)
set(EXTENSION_HOMEPAGE "http://www.carmacenter.org")
set(EXTENSION_CATEGORY "CARMA")
set(EXTENSION_STATUS "Beta")
set(EXTENSION_ICONURL "http://www.na-mic.org/Wiki/images/a/ac/Carma-logo-128.png")
set(EXTENSION_CONTRIBUTORS "Alan Morris (CARMA), Danny Perry (CARMA), Josh Cates (CARMA), Rob MacLeod (CARMA), Eugene Kholmovski (CARMA)")
set(EXTENSION_DESCRIPTION "The CARMA extension includes modules for the Utah Afib DBP")
set(EXTENSION_SCREENSHOTURLS "http://wiki.na-mic.org/Wiki/images/d/d1/Carma_afib_inhomogeneity_correction.png http://wiki.na-mic.org/Wiki/images/7/7a/Carma_afib_auto_scar.png")

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()


set(MODULE_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  )

#-----------------------------------------------------------------------------
# Inhomogeneity Correction Module
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
set(MODULE_NAME "CarmaInhomogeneityCorrection")

#-----------------------------------------------------------------------------
slicerMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/NAMICLogo.h
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  )


#-----------------------------------------------------------------------------
# Automatic Left Atrial Scar Module
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
set(MODULE_NAME CarmaAutomaticLeftAtrialScar)

#-----------------------------------------------------------------------------
slicerMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/NAMICLogo.h
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()

