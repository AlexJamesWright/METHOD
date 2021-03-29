project(METHOD)
add_library(Cmin "")

# HDF5 CMake scripts don't work properly. This is hard-coded.
# Commented out for Travis
# include_directories(/usr/include/hdf5/openmpi)

find_package(MPI)
if(MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
endif()

target_sources(
    Cmin
    PRIVATE
      ${CMAKE_CURRENT_LIST_DIR}/Src/TimeIntegrators/backRKArgs.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/TimeIntegrators/backwardsRK.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/BoundaryConditions/boundaryConds.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/Flux/fluxVectorSplitting.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/TimeIntegrators/IMEX2Args.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/TimeIntegrators/IMEX3Args.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/InitialData/initFunc.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/TimeIntegrators/RK2.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/TimeIntegrators/rkSplit.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/DataIO/serialSaveData.cc
      # ${CMAKE_CURRENT_LIST_DIR}/Src/DataIO/serialSaveDataHDF5.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/Simulation/serialEnv.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/Simulation/simData.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/Simulation/simulation.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/Physics/srmhd.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/Physics/srrmhd.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/Physics/REGIME.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/TimeIntegrators/SSP2.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/TimeIntegrators/SSP3.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/TimeIntegrators/SSP2322.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/Physics/twoFluidEMHD.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/Flux/weno.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/Flux/wenoUpwinds.cc
      ${CMAKE_CURRENT_LIST_DIR}/Src/TimeIntegrators/RKPlus.cc

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
    ${CMAKE_CURRENT_LIST_DIR}/Include/TimeIntegrators/backRKArgs.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/TimeIntegrators/backwardsRK.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/BoundaryConditions/boundaryConds.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/Flux/flux.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/Flux/fluxVectorSplitting.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/TimeIntegrators/IMEX2Args.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/TimeIntegrators/IMEX3Args.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/InitialData/initFunc.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/Physics/model.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/TimeIntegrators/RK2.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/TimeIntegrators/rkSplit.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/DataIO/saveData.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/DataIO/serialSaveData.h
    # ${CMAKE_CURRENT_LIST_DIR}/Include/DataIO/serialSaveDataHDF5.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/Simulation/platformEnv.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/Simulation/serialEnv.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/Simulation/simData.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/Simulation/simulation.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/Physics/srmhd.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/Physics/srrmhd.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/TimeIntegrators/SSP2.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/TimeIntegrators/SSP3.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/TimeIntegrators/SSP2322.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/Physics/modelExtension.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/Physics/REGIME.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/TimeIntegrators/timeInt.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/Physics/twoFluidEMHD.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/Flux/wenoUpwinds.h
    ${CMAKE_CURRENT_LIST_DIR}/Include/TimeIntegrators/RKPlus.h

    ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Include/cminpack.h
    ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Include/cminpackP.h
    ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Include/minpack.h
    ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Include/minpackP.h
    )

target_include_directories(
    Cmin
    PUBLIC
      ${CMAKE_CURRENT_LIST_DIR}/Include/Simulation
      ${CMAKE_CURRENT_LIST_DIR}/Include/Flux      
      ${CMAKE_CURRENT_LIST_DIR}/Include/TimeIntegrators
      ${CMAKE_CURRENT_LIST_DIR}/Include/DataIO
      ${CMAKE_CURRENT_LIST_DIR}/Include/Physics
      ${CMAKE_CURRENT_LIST_DIR}/Include/InitialData
      ${CMAKE_CURRENT_LIST_DIR}/Include/BoundaryConditions
      ${CMAKE_CURRENT_LIST_DIR}/CminpackLibrary/Include
    )
