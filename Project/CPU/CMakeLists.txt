project(METHOD)
add_library(Cmin "" ../../Examples/h5test.cpp)

# Added HDF5 to Cmake for Clion
find_package(HDF5)
if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIR})
    set(_hdf5_libs hdf5 hdf5_cpp)
endif()

find_package(MPI)
if(MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
endif()

target_sources(
    Cmin
    PRIVATE
      ${CMAKE_CURRENT_LIST_DIR}/Src/backRKArgs.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/backwardsRK.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/boundaryConds.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/fluxVectorSplitting.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/IMEX2Args.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/IMEX3Args.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/initFunc.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/RK2.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/rkSplit.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/serialSaveData.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/serialSaveDataHDF5.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/serialEnv.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/simData.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/simulation.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/srmhd.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/srrmhd.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/REGIME.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/SSP2.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/SSP3.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/SSP2322.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/twoFluidEMHD.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/weno.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/wenoUpwinds.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/RKPlus.cc

      ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Src/dogleg.cc
      ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Src/dpmpar.cc
      ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Src/enorm.cc
      ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Src/fdjac1.cc
      ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Src/hybrd1.cc
      ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Src/hybrd.cc
      ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Src/qform.cc
      ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Src/qrfac.cc
      ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Src/r1mpyq.cc
      ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Src/r1updt.cc
    PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/Include/backRKArgs.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/backwardsRK.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/boundaryConds.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/flux.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/fluxVectorSplitting.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/IMEX2Args.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/IMEX3Args.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/initFunc.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/model.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/RK2.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/rkSplit.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/saveData.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/serialSaveData.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/serialSaveDataHDF5.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/platformEnv.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/serialEnv.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/simData.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/simulation.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/srmhd.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/srrmhd.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/SSP2.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/SSP3.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/SSP2322.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/modelExtension.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/REGIME.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/timeInt.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/twoFluidEMHD.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/wenoUpwinds.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/RKPlus.h

    ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Include/cminpack.h
    ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Include/cminpackP.h
    ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Include/minpack.h
    ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Include/minpackP.h
    )

target_include_directories(
    Cmin
    PUBLIC
      ${CMAKE_CURRENT_LIST_DIR}/Include
      ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Include
    )
