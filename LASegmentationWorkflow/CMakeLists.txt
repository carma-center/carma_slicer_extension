cmake_minimum_required(VERSION 2.8.7)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# LA Segmentation Workflow
#-----------------------------------------------------------------------------
set(MODULE_NAME LASegmentationWorkflow)
#-----------------------------------------------------------------------------

set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  )

set(MODULE_PYTHON_RESOURCES
  LASegmentationResources/${MODULE_NAME}SelectDataStep.py
  LASegmentationResources/${MODULE_NAME}EndoSegmentationStep.py
  LASegmentationResources/${MODULE_NAME}ImageRegistrationStep.py
  #LASegmentationResources/${MODULE_NAME}EndoCleanupStep.py
  LASegmentationResources/${MODULE_NAME}Step.py
  LASegmentationResources/__init__.py
  LASegmentationResources/Helper.py
  )

slicerMacroBuildScriptedModule(
  NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
  #DESTINATION_DIR ${CMAKE_BINARY_DIR}/${Slicer_QTSCRIPTEDMODULES_LIB_DIR}
  #INSTALL_DIR ${Slicer_INSTALL_QTSCRIPTEDMODULES_LIB_DIR}
  )
