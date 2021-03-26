#include "gtest/gtest.h"
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srrmhd.h"
#include "boundaryConds.h"
#include "SSP2.h"
#include "SSP3.h"
#include "saveData.h"
#include "fluxVectorSplitting.h"
#include "serialSaveData.h"
#include <cstdlib>
#include <cstdio>

TEST(SSP2, IMEX2ConsistentWithSerialVersion)
{
  /*
    Run a resistive MHD problem to test the IMEX2 scheme. This test will run
    the simulation and save the output data in the TestData directory, ready
    to be compared to the serial output.
  */
  double sigma(0);

  SerialEnv env(0, NULL, 1, 1, 1);
  Data data(64, 16, 0, 0, 1, 0, 1, 0, 1, 0.05, &env,
            0.5, 4, 4.0/3.0, sigma);

  // Choose particulars of simulation
  SRRMHD model(&data);
  FVS fluxMethod(&data, &model);
  Outflow bcs(&data);
  Simulation sim(&data, &env);
  BrioWuSingleFluid init(&data);
  SSP2 timeInt(&data, &model, &bcs, &fluxMethod);
  SerialSaveData save(&data, &env);

  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  sim.evolve();

  // Save data in test directory
  strcpy(save.dir, "../TestData/GPU");
  strcpy(save.app, "SSP2");
  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}

TEST(SSP2FlowKHSingleFluid, IMEX2ConsistentWithSerialVersion)
{
  /*
    Run a resistive MHD problem to test the IMEX2 scheme. This test will run
    the simulation and save the output data in the TestData directory, ready
    to be compared to the serial output.
  */
  double sigma(0);

  SerialEnv env(0, NULL, 1, 1, 1);
  Data data(64, 16, 0, -0.5, 0.5, -1, 1, 0, 1, 0.05, &env,
            0.5, 4, 4.0/3.0, sigma);

  // Choose particulars of simulation
  SRRMHD model(&data);
  FVS fluxMethod(&data, &model);
  Flow bcs(&data);
  Simulation sim(&data, &env);
  KHInstabilitySingleFluid init(&data);
  SSP2 timeInt(&data, &model, &bcs, &fluxMethod);
  SerialSaveData save(&data, &env);

  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  sim.evolve();

  // Save data in test directory
  strcpy(save.dir, "../TestData/GPU");
  strcpy(save.app, "SSP2FlowKHSingleFluid");
  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}

#if 0
TEST(SSP2FlowKHSingleFluid, IMEX2ConsistentWithSerialVersion)
{
  /*
    Run a resistive MHD problem to test the IMEX2 scheme. This test will run
    the simulation and save the output data in the TestData directory, ready
    to be compared to the serial output.
  */
  const double MU(1000); 
  double sigma(300);
  double cp(1.0);
  double mu1(-MU);
  double mu2(MU); 
  int nx(256);
  int ny(512); 
  double endTime(0.01);

  Data data(nx, ny, 0, -0.5, 0.5, -1, 1, 0, 1, endTime,
            0.1, 4, 4.0/3.0, sigma, cp, mu1, mu2);

  // Choose particulars of simulation
  SRRMHD model(&data);
  FVS fluxMethod(&data, &model);
  Simulation sim(&data);
  KHInstabilitySingleFluid init(&data);
  Flow bcs(&data);
  SSP2 timeInt(&data, &model, &bcs, &fluxMethod);
  SaveData save(&data);

  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  sim.evolve();

  // Save data in test directory
  strcpy(save.dir, "../TestData/GPU");
  strcpy(save.app, "SSP2FlowKHSingleFluid");
  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}
#endif

TEST(SSP3, IMEX3ConsistentWithSerialVersion)
{
  /*
    Run a resistive MHD problem to test the IMEX3 scheme. This test will run
    the simulation and save the output data in the TestData directory, ready
    to be compared with the parallel output.
  */
  double sigma(0);

  SerialEnv env(0, NULL, 1, 1, 1);
  Data data(64, 16, 0, 0, 1, 0, 1, 0, 1, 0.05, &env,
            0.5, 4, 4.0/3.0, sigma);

  // Choose particulars of simulation
  SRRMHD model(&data);
  FVS fluxMethod(&data, &model);
  Outflow bcs(&data);
  Simulation sim(&data, &env);
  BrioWuSingleFluid init(&data);
  SSP3 timeInt(&data, &model, &bcs, &fluxMethod);
  SerialSaveData save(&data, &env);

  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  sim.evolve();

  // Save data in test directory
  strcpy(save.dir, "../TestData/GPU");
  strcpy(save.app, "SSP3");
  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}
