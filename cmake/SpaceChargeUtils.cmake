macro(project_add_gtest TESTNAME FILES LIBRARIES TEST_WORKING_DIRECTORY)
  add_executable(${TESTNAME} ${FILES})
  target_link_libraries(${TESTNAME} gtest gmock gtest_main ${LIBRARIES})
  gtest_discover_tests(
    ${TESTNAME}
    WORKING_DIRECTORY ${TEST_WORKING_DIRECTORY}
    PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${TEST_WORKING_DIRECTORY}")
  set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()