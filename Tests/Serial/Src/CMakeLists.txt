add_executable(
    unit_tests
    test_boundaryConds.cc
    test_cminpack.cc
    test_initFunc.cc
    # test_resistiveSGM.cc
    test_simulation.cc
    test_srmhd.cc
    test_srrmhd.cc
    test_twoFluidEMHD.cc
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
