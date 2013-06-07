cmake_minimum_required(VERSION 2.8.7)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Cardiac MRI Toolkit Wizard
#-----------------------------------------------------------------------------
set(MODULE_NAME CMRToolkitWizard)
#-----------------------------------------------------------------------------

set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  )

set(MODULE_PYTHON_RESOURCES
  Resources/${MODULE_NAME}SelectDataStep.py
  Resources/${MODULE_NAME}EndoSegmentationStep.py
  Resources/${MODULE_NAME}AxialDilateStep.py
  Resources/${MODULE_NAME}BooleanRemoveStep.py
  Resources/${MODULE_NAME}WallCleanupStep.py
  Resources/${MODULE_NAME}AntrumCutStep.py
  Resources/${MODULE_NAME}EndoNVIsosurfaceStep.py
  Resources/${MODULE_NAME}AutoScarStep.py
  Resources/${MODULE_NAME}Step.py
  Resources/__init__.py
  Resources/Helper.py
  )

__slicerMacroBuildScriptedModule__(
  TARGET_NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
  DESTINATION_DIR ${CMAKE_BINARY_DIR}/${Slicer_QTSCRIPTEDMODULES_LIB_DIR}
  INSTALL_DIR ${Slicer_INSTALL_QTSCRIPTEDMODULES_LIB_DIR}
  )
