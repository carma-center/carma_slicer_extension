
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
 # PVAtrumCut
 ##############################################################################
 set(CLP "PVAntrumCut")

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
   --Endo_Layer ${INPUT}/shrinkwrap-endo.nrrd --Wall_Layer ${INPUT}/shrinkwrap-wall.nrrd --Endo_Layer_No_Veins ${TEMP}/${CLP}Test.nrrd
   )
 set_property(TEST ${testname} PROPERTY LABELS ${CLP})
