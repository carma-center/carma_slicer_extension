 #-----------------------------------------------------------------------------
 set(Launcher_Command ${Slicer_LAUNCH_COMMAND})
 set(BASELINE ${CMAKE_CURRENT_SOURCE_DIR}/../../Data/Baseline)
 set(INPUT ${CMAKE_CURRENT_SOURCE_DIR}/../../Data/Input)
 set(TEMP "${CMAKE_BINARY_DIR}/Testing/Temporary")


 ##############################################################################
 # CarmaAutomaticLeftAtrialScar
 ##############################################################################
 set(CLP "CarmaAutomaticLeftAtrialScar")

 #-----------------------------------------------------------------------------
 add_executable(${CLP}Test ${CLP}Test.cxx)
 target_link_libraries(${CLP}Test ${CLP}Lib)
 set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})

 #-----------------------------------------------------------------------------
 set(testname ${CLP}Test)
 add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
   --compare ${BASELINE}/${CLP}Test.nrrd
   ${TEMP}/${CLP}Test.nrrd
   ModuleEntryPoint
   ${INPUT}/autoscar-data.nrrd ${INPUT}/autoscar-wall.nrrd ${INPUT}/autoscar-endo.nrrd ${TEMP}/${CLP}Test.nrrd
   )
 set_property(TEST ${testname} PROPERTY LABELS ${CLP})

 ##############################################################################
 # CarmaInhomogeneityCorrection
 ##############################################################################
 set(CLP "CarmaInhomogeneityCorrection")

 #-----------------------------------------------------------------------------
 add_executable(${CLP}Test ${CLP}Test.cxx)
 target_link_libraries(${CLP}Test ${CLP}Lib)
 set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})

 #-----------------------------------------------------------------------------
 set(testname ${CLP}Test)
 add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
   --compare ${BASELINE}/${CLP}Test.nrrd
   ${TEMP}/${CLP}Test.nrrd
   ModuleEntryPoint
   ${INPUT}/inhomogeneity-data.nrrd ${INPUT}/inhomogeneity-endo.nrrd ${INPUT}/inhomogeneity-roi.nrrd ${TEMP}/${CLP}Test.nrrd
   )
 set_property(TEST ${testname} PROPERTY LABELS ${CLP})

 ##############################################################################
 # CarmaAxialDilate
 ##############################################################################
 set(CLP "CarmaAxialDilate")

 #-----------------------------------------------------------------------------
 add_executable(${CLP}Test ${CLP}Test.cxx)
 target_link_libraries(${CLP}Test ${CLP}Lib)
 set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})

 #-----------------------------------------------------------------------------
 set(testname ${CLP}Test)
 add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
   --compare ${BASELINE}/${CLP}Test.nrrd
   ${TEMP}/${CLP}Test.nrrd
   ModuleEntryPoint
   ${INPUT}/dilate-endo.nrrd ${TEMP}/${CLP}Test.nrrd
   )
 set_property(TEST ${testname} PROPERTY LABELS ${CLP})

 ##############################################################################
 # CarmaPVAntrumCut
 ##############################################################################
 set(CLP "CarmaPVAntrumCut")

 #-----------------------------------------------------------------------------
 add_executable(${CLP}Test ${CLP}Test.cxx)
 target_link_libraries(${CLP}Test ${CLP}Lib)
 set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})

 #-----------------------------------------------------------------------------
 set(testname ${CLP}Test)
 add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
   --compare ${BASELINE}/${CLP}Test.nrrd
   ${TEMP}/${CLP}Test.nrrd
   ModuleEntryPoint
   ${INPUT}/shrinkwrap-endo.nrrd ${INPUT}/shrinkwrap-wall.nrrd ${TEMP}/${CLP}Test.nrrd
   )
 set_property(TEST ${testname} PROPERTY LABELS ${CLP})

 ##############################################################################
 # CarmaIsolatedConnectedFilter
 ##############################################################################
 set(CLP "CarmaIsolatedConnectedFilter")

 #-----------------------------------------------------------------------------
 add_executable(${CLP}Test ${CLP}Test.cxx)
 target_link_libraries(${CLP}Test ${CLP}Lib)
 set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})

 #-----------------------------------------------------------------------------
 set(testname ${CLP}Test)
 add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
   ModuleEntryPoint
   --inputImage ${INPUT}/segmentation-data.nrrd 
   --outputImage ${TEMP}/${CLP}Test.nrrd 
   --lowerThreshold 150 
   --seeds1 4.42729,-6.17517,34.6671 
   --seeds2 72.593,21.5171,32.1671
   )
 set_property(TEST ${testname} PROPERTY LABELS ${CLP})

 ##############################################################################
 # CarmaConnectedThresholdFilter
 ##############################################################################
 set(CLP "CarmaConnectedThresholdFilter")

 #-----------------------------------------------------------------------------
 add_executable(${CLP}Test ${CLP}Test.cxx)
 target_link_libraries(${CLP}Test ${CLP}Lib)
 set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})

 #-----------------------------------------------------------------------------
 set(testname ${CLP}Test)
 add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
   ModuleEntryPoint
   --inputImage ${INPUT}/segmentation-data.nrrd 
   --outputImage ${TEMP}/${CLP}Test.nrrd 
   --lowerThreshold 150 
   --upperThreshold 310 
   --seeds 16.99,-0.290722,27.1671 
   --seeds -28.8485,-3.64476,27.1671
   )
 set_property(TEST ${testname} PROPERTY LABELS ${CLP})

##############################################################################
 # CarmaBooleanREMOVE
 ##############################################################################
 set(CLP "CarmaBooleanREMOVEFilter")

 #-----------------------------------------------------------------------------
 add_executable(${CLP}Test ${CLP}Test.cxx)
 target_link_libraries(${CLP}Test ${CLP}Lib)
 set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})

 #-----------------------------------------------------------------------------
 set(testname ${CLP}Test)
 add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
   ModuleEntryPoint
   ${INPUT}/boolean-endodilate.nrrd ${INPUT}/boolean-endo.nrrd ${TEMP}/${CLP}Test.nrrd
   )
 set_property(TEST ${testname} PROPERTY LABELS ${CLP})
