cmake_minimum_required(VERSION 2.8.7)

set(EXTENSION_NAME carma)
set(EXTENSION_HOMEPAGE "http://www.carmacenter.org")
set(EXTENSION_CATEGORY "CARMA")
set(EXTENSION_STATUS "Beta")
set(EXTENSION_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/LoadableExtensionTemplate/Resources/Icons/LoadableExtensionTemplate.png?revision=19437&view=co")
set(EXTENSION_CONTRIBUTORS "Alan Morris (CARMA), Danny Perry (CARMA), Josh Cates (CARMA), Rob MacLeod (CARMA), Eugene Kholmovski (CARMA)")
set(EXTENSION_DESCRIPTION "The CARMA extension includes modules for the Utah Afib DBP")
set(EXTENSION_SCREENSHOTURLS "http://wiki.na-mic.org/Wiki/images/d/d1/Carma_afib_inhomogeneity_correction.png http://wiki.na-mic.org/Wiki/images/7/7a/Carma_afib_auto_scar.png")

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

set(cli_module_target_libraries
  ${ITK_LIBRARIES}
)



#-----------------------------------------------------------------------------
# Inhomogeneity Correction Module
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
set(MODULE_NAME "Inhomogeneity Correction") # Do not use 'project()'
set(MODULE_FILE_NAME "CarmaInhomogeneityCorrection")

#-----------------------------------------------------------------------------
set(MODULE_HOMEPAGE "http://wiki.na-mic.org/Wiki/index.php/DBP3:Utah:SlicerModuleInhomogeneity")
set(MODULE_CATEGORY "CARMA")
set(MODULE_STATUS "Beta")
set(MODULE_CONTRIBUTORS "Alan Morris (CARMA), Eugene Kholmovski (CARMA), Josh Cates (CARMA), Danny Perry (CARMA), Rob MacLeod (CARMA)")
set(MODULE_DESCRIPTION "Inhomogeneity Correction for the CARMA Afib project")
set(MODULE_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/carma_inhomogeneity_correction/carma_inhomogeneity_correction.png?revision=19437&view=co")
set(MODULE_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/a/ab/Slicer-r19441-inhomogeneity_correction-screenshot.png http://wiki.slicer.org/slicerWiki/images/1/1e/Slicer-r19441-inhomogeneity_correction-screenshot-2.png")
set(MODULE_ACKNOWLEDGEMENTS "This work was supported by ...")
set(MODULE_MAJOR_VERSION 0)
set(MODULE_MINOR_VERSION 1)
set(MODULE_PATCH_VERSION 0)



#-----------------------------------------------------------------------------
# Configure XML description
configure_file(
  ${MODULE_FILE_NAME}.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_FILE_NAME}.xml
  @ONLY
  )


#-----------------------------------------------------------------------------
slicerMacroBuildCLI(
  NAME ${MODULE_FILE_NAME}
  #LOGO_HEADER ${Slicer_SOURCE_DIR}/Resources/NAMICLogo.h
  TARGET_LIBRARIES ${cli_module_target_libraries}
  INCLUDE_DIRECTORIES ${cli_module_include_directories}
  ADDITIONAL_SRCS ${cli_module_SRCS}
  #EXECUTABLE_ONLY
  )


#-----------------------------------------------------------------------------
# Automatic Left Atrial Scar Module
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
set(MODULE_NAME "Automatic Left Atrial Scar") # Do not use 'project()'
set(MODULE_FILE_NAME "CarmaAutomaticLeftAtrialScar")

#-----------------------------------------------------------------------------
set(MODULE_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/Snapshot/Developers/Tutorials/ExtensionWriting")
set(MODULE_CATEGORY "CARMA")
set(MODULE_STATUS "Beta")
set(MODULE_CONTRIBUTORS "Danny Perry (CARMA), Alan Morris (CARMA), Josh Cates (CARMA), Rob MacLeod (CARMA)")
set(MODULE_DESCRIPTION "Determine Left Atrial Scar Automatically given the wall segmentation, in a post-ablation LGE-MRI. For details see: Perry, D. et al. "Automatic classification of scar tissue...", Proceedings of SPIE Medical Imaging: Computer Aided Diagnosis, Feb 2012.")
set(MODULE_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/CarmaAutomaticLeftAtrialScar/CarmaAutomaticLeftAtrialScar.png?revision=19437&view=co")
set(MODULE_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/a/ab/Slicer-r19441-CarmaAutomaticLeftAtrialScar-screenshot.png http://wiki.slicer.org/slicerWiki/images/1/1e/Slicer-r19441-CarmaAutomaticLeftAtrialScar-screenshot-2.png")
set(MODULE_ACKNOWLEDGEMENTS "This work was supported by ...")
set(MODULE_MAJOR_VERSION 0)
set(MODULE_MINOR_VERSION 1)
set(MODULE_PATCH_VERSION 0)

#-----------------------------------------------------------------------------
# Configure XML description
configure_file(
  ${MODULE_FILE_NAME}.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_FILE_NAME}.xml
  @ONLY
  )


#-----------------------------------------------------------------------------
slicerMacroBuildCLI(
  NAME ${MODULE_FILE_NAME}
  #LOGO_HEADER ${Slicer_SOURCE_DIR}/Resources/NAMICLogo.h
  TARGET_LIBRARIES ${cli_module_target_libraries}
  INCLUDE_DIRECTORIES ${cli_module_include_directories}
  ADDITIONAL_SRCS ${cli_module_SRCS}
  #EXECUTABLE_ONLY
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()
