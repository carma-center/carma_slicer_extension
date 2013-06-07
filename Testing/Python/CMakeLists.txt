
#-----------------------------------------------------------------------------
#set(KIT_UNITTEST_SCRIPTS)
#SlicerMacroConfigureGenericPythonModuleTests("CarmaRegistration" KIT_UNITTEST_SCRIPTS)

#-----------------------------------------------------------------------------
#foreach(script_name ${KIT_UNITTEST_SCRIPTS})
#  slicer_add_python_unittest(
#    SCRIPT ${script_name}
#    SLICER_ARGS --no-main-window --disable-cli-modules --disable-loadable-modules
#                --additional-module-path ${CMAKE_BINARY_DIR}/${Slicer_QTSCRIPTEDMODULES_LIB_DIR}
#    TESTNAME_PREFIX nomainwindow_
#    )
#endforeach()

#if(Slicer_USE_QtTesting AND Slicer_USE_PYTHONQT)

# slicer_add_python_unittest(SCRIPT CardiacMRIToolkitSelfTest.py)

# set(KIT_PYTHON_SCRIPTS
#    CardiacMRIToolkitSelfTest.py
#    )

# set(KIT_PYTHON_RESOURCES
#    Resources/UI/ScenePerformance.ui
    )

# __slicerMacroBuildScriptedModule__(
#    TARGET_NAME ApplicationSelfTests
#    SCRIPTS "${KIT_PYTHON_SCRIPTS}"
#    RESOURCES "${KIT_PYTHON_RESOURCES}"
#    DESTINATION_DIR ${CMAKE_BINARY_DIR}/${Slicer_QTSCRIPTEDMODULES_LIB_DIR}
#    INSTALL_DIR ${Slicer_INSTALL_QTSCRIPTEDMODULES_LIB_DIR}
#    NO_INSTALL_SUBDIR
#    )

#endif()