include(Catch2)

macro (nwgraph_add_hpx_exe EXENAME)
  add_executable(${EXENAME}.exe ${EXENAME}.cpp)
  target_link_libraries(${EXENAME}.exe nwgraph_hpx)
endmacro (nwgraph_add_hpx_exe)

macro (nwgraph_add_hpx_test TESTNAME)
  add_executable(${TESTNAME}.exe ${TESTNAME}.cpp ../hpx_main.cpp)
  target_link_libraries(${TESTNAME}.exe nwgraph nwgraph_hpx Catch2::Catch2)
  target_compile_definitions(${TESTNAME}.exe PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data/")
  target_compile_definitions(${TESTNAME}.exe PRIVATE DATA_FILE="${CMAKE_SOURCE_DIR}/data/karate.mtx")
  add_test(NAME ${TESTNAME}_test COMMAND ${TESTNAME}.exe)
  #  catch_discover_tests(${TESTNAME})
endmacro (nwgraph_add_hpx_test)

nwgraph_add_hpx_test(partitioned_adj_io_test)
