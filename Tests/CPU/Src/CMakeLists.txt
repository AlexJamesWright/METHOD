add_executable(
    unit_tests
    test_twoFluidEMHD.cc
    test_srmhd.cc
    test_srrmhd.cc
    test_simulation.cc
    test_initFunc.cc
    test_cminpack.cc
    test_boundaryConds.cc
    test_REGIME.cc
    )

target_link_libraries(
    unit_tests
    gtest_main
    Cmin
    )

add_test(
  NAME
    unit
  COMMAND
    ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/unit_tests
  )
