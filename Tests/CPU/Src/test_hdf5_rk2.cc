#include "gtest/gtest.h"
#include "srrmhd.h"
#include "srmhd.h"
#include "simulation.h"
#include "simData.h"
#include "saveData.h"
#include "serialSaveDataHDF5.h"
#include "initFunc.h"
#include "RK2.h"
#include "fluxVectorSplitting.h"
#include "serialEnv.h"
#include <cstdlib>

#if 1
// RKOTVSingleFluidPeriodic
TEST(RK2OutputConsistentWithSerial, RK2SrmhdPeriodicOTVSF)
{

  /*
  The following was used to gather data to compare the parallel
  version with. No tests are run in the serial version of this test
  */

  double cfl(0.6);
  int Ng(4);
  double gamma(2.0);

  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env, cfl, Ng, gamma);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  Periodic bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveDataHDF5 save(&d, &env, "../TestData/SerialHDF5/RK2SrmhdPeriodicOTVSF", SerialSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // sim.evolve();
  sim.updateTime();
  // sim.updateTime();


  // Save data in test directory
  // This currently needs to be set in the save() function above as well
  strcpy(save.dir, "../TestData/SerialHDF5");
  strcpy(save.app, "RK2SrmhdPeriodicOTVSF");

  save.saveAll();
}

// RKOTVSingleFluidOutflow
TEST(RK2OutputConsistentWithSerial, RK2SrmhdOutflowOTVSF)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  double cfl(0.6);
  int Ng(4);
  double gamma(2.0);

  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env, cfl, Ng, gamma);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  Outflow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveDataHDF5 save(&d, &env, "../TestData/SerialHDF5/RK2SrmhdOutflowOTVSF", SerialSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // sim.evolve();
  sim.updateTime();
  // sim.updateTime();


  // Save data in test directory
  strcpy(save.dir, "../TestData/SerialHDF5");
  strcpy(save.app, "RK2SrmhdOutflowOTVSF");

  save.saveAll();
}


// RKRandomInstabilitySingleFluid
TEST(RK2OutputConsistentWithSerial, RK2SrmhdOutflowKHRandomInstabilitySF)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  int seed(10);
  double cfl(0.6);
  int Ng(4);
  double gamma(4.0/3.0);
  double sigma(10);
  double cp(1.0);
  double mu1(-100);
  double mu2(100);
  int frameSkip(10);

  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env, cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  Outflow bcs(&d);
  Simulation sim(&d, &env);
  KHRandomInstabilitySingleFluid init(&d, 1, seed);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveDataHDF5 save(&d, &env, "../TestData/SerialHDF5/RK2SrmhdOutflowKHRandomInstabilitySF", SerialSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/SerialHDF5");
  strcpy(save.app, "RK2SrmhdOutflowKHRandomInstabilitySF");

  save.saveAll();
}

TEST(RK2OutputConsistentWithSerial, RK2SrmhdPeriodicKHRandomInstabilitySF)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  int seed(10);
  double cfl(0.6);
  int Ng(4);
  double gamma(4.0/3.0);
  double sigma(10);
  double cp(1.0);
  double mu1(-100);
  double mu2(100);
  int frameSkip(10);

  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env, cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip);
  SRMHD model(&d);
  Weno3 weno(&d);
  FVS fluxMethod(&d, &weno, &model);
  Periodic bcs(&d);
  Simulation sim(&d, &env);
  KHRandomInstabilitySingleFluid init(&d, 1, seed);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveDataHDF5 save(&d, &env, "../TestData/SerialHDF5/RK2SrmhdPeriodicKHRandomInstabilitySF", SerialSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/SerialHDF5");
  strcpy(save.app, "RK2SrmhdPeriodicKHRandomInstabilitySF");

  save.saveAll();
}

TEST(RK2OutputConsistentWithSerial, RK2SrmhdFlowKHRandomInstabilitySF)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  int seed(10);
  double cfl(0.6);
  int Ng(4);
  double gamma(4.0/3.0);
  double sigma(10);
  double cp(1.0);
  double mu1(-100);
  double mu2(100);
  int frameSkip(10);

  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env, cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  KHRandomInstabilitySingleFluid init(&d, 1, seed);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveDataHDF5 save(&d, &env, "../TestData/SerialHDF5/RK2SrmhdFlowKHRandomInstabilitySF", SerialSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/SerialHDF5");
  strcpy(save.app, "RK2SrmhdFlowKHRandomInstabilitySF");

  save.saveAll();
}

// BrioWuSingleFluid

TEST(RK2OutputConsistentWithSerial, RK2SrmhdOutflowBrioWuSF)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  Outflow bcs(&d);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveDataHDF5 save(&d, &env, "../TestData/SerialHDF5/RK2SrmhdOutflowBrioWuSF", SerialSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/SerialHDF5");
  strcpy(save.app, "RK2SrmhdOutflowBrioWuSF");

  save.saveAll();
}

TEST(RK2OutputConsistentWithSerial, RK2SrmhdPeriodicBrioWuSF)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  Periodic bcs(&d);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveDataHDF5 save(&d, &env, "../TestData/SerialHDF5/RK2SrmhdPeriodicBrioWuSF", SerialSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/SerialHDF5");
  strcpy(save.app, "RK2SrmhdPeriodicBrioWuSF");

  save.saveAll();
}

TEST(RK2OutputConsistentWithSerial, RK2SrmhdFlowBrioWuSF)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveDataHDF5 save(&d, &env, "../TestData/SerialHDF5/RK2SrmhdFlowBrioWuSF", SerialSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/SerialHDF5");
  strcpy(save.app, "RK2SrmhdFlowBrioWuSF");

  save.saveAll();
}
#endif

#if 0
TEST(RK2OutputConsistentWithSerial, RK2SrrmhdOutflowOTVortexSingleFluidRK2)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  SRRMHD model(&d);
  FVS fluxMethod(&d, &model);
  Outflow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveDataHDF5 save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/SerialHDF5");
  strcpy(save.app, "RK2SrrmhdOutflowOTVortexSingleFluid");

  save.saveAll();
}
#endif
