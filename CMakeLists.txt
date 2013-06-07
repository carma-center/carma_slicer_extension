cmake_minimum_required(VERSION 2.8.7)

#-----------------------------------------------------------------------------
set(EXTENSION_NAME "Cardiac_MRI_Toolkit")
set(EXTENSION_HOMEPAGE "http://capulet.med.utah.edu/namic")
set(EXTENSION_CATEGORY "Cardiac MRI toolkit")
set(EXTENSION_STATUS "Beta")
set(EXTENSION_ICONURL "http://capulet.med.utah.edu/namic/images/CMRTK-logo.png")
set(EXTENSION_CONTRIBUTORS "Alan Morris (CARMA), Salma Bengali (CARMA), Josh Cates (CARMA), Rob MacLeod (CARMA), Eugene Kholmovski (CARMA)")
set(EXTENSION_DESCRIPTION "The Cardiac MRI Toolkit includes modules for the Utah Afib DBP")
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
set(MODULE_NAME "CMRToolkitInhomogeneityCorrection")

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/NamicCarmaLogo.h
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  )

#-----------------------------------------------------------------------------
# Axial Dilate Module
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
set(MODULE_NAME "CMRToolkitAxialDilate")

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/NamicCarmaLogo.h
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  )

#-----------------------------------------------------------------------------
# Registration Module
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
set(MODULE_NAME "CMRToolkitRegistration")

#-----------------------------------------------------------------------------
slicerMacroBuildScriptedModule(
  NAME CMRToolkitRegistration
  SCRIPTS "CMRToolkitRegistration.py"
  RESOURCES ""
  )


#-----------------------------------------------------------------------------
# Automatic Left Atrial Scar Module
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
set(MODULE_NAME CMRToolkitAutomaticLeftAtrialScar)

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/NamicCarmaLogo.h
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  )

#-----------------------------------------------------------------------------
# Registration Module using BRAINSFit
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
set(MODULE_NAME "CMRToolkitRegistrationBRAINSFit")

#-----------------------------------------------------------------------------
slicerMacroBuildScriptedModule(
  NAME CMRToolkitRegistrationBRAINSFit
  SCRIPTS "CMRToolkitRegistrationBRAINSFit.py"
  RESOURCES ""
  #DESTINATION_DIR ${CMAKE_BINARY_DIR}/${Slicer_QTSCRIPTEDMODULES_LIB_DIR}
  #INSTALL_DIR ${Slicer_INSTALL_QTSCRIPTEDMODULES_LIB_DIR}
  #NO_INSTALL_SUBDIR
  )

#-----------------------------------------------------------------------------
# PV Antrum Cut Module
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
set(MODULE_NAME CMRToolkitPVAntrumCut)

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/NamicCarmaLogo.h
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  )

#-----------------------------------------------------------------------------
# Connected Threshold Filter Module
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
set(MODULE_NAME CMRToolkitConnectedThresholdFilter)

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/NamicCarmaLogo.h
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  )

#-----------------------------------------------------------------------------
# Isolated Connected Filter Module
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
set(MODULE_NAME CMRToolkitIsolatedConnectedFilter)

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/NamicCarmaLogo.h
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  )

#-----------------------------------------------------------------------------
# Boolean REMOVE Filter Module
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
set(MODULE_NAME CMRToolkitBooleanREMOVEFilter)

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/NamicCarmaLogo.h
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  )

#-----------------------------------------------------------------------------
# Add LA Segmentation Using Graph Cuts Module (In progress)
#-----------------------------------------------------------------------------
# add_subdirectory(CMRToolkitLASegmentationGraphCuts)

#-----------------------------------------------------------------------------
# Add Cardiac MRI Toolkit Wizard
#-----------------------------------------------------------------------------
add_subdirectory(CMRToolkitWizard)

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()
