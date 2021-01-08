#include "gtest/gtest.h"
#include "srrmhd.h"
#include "srmhd.h"
#include "simulation.h"
#include "simData.h"
#include "parallelSaveDataHDF5.h"
#include "parallelBoundaryConds.h"
#include "initFunc.h"
#include "RK2.h"
#include "fluxVectorSplitting.h"
#include "parallelEnv.h"
#include <cstdlib>

/*
 Assumptions:
   RKRandomInstabilitySingleFluid is tested in 2D only
   BrioWuSingleFluid is tested in 1D only
*/


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

  ParallelEnv env(0, NULL, 2, 2, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env, cfl, Ng, gamma);
  SRMHD model(&d);
  Weno3 weno(&d);
  FVS fluxMethod(&d, &weno, &model);
  ParallelPeriodic bcs(&d, &env);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveDataHDF5 save(&d, &env, "../TestData/CPUHDF5/RK2SrmhdPeriodicOTVSF", ParallelSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // sim.evolve();
  sim.updateTime();
  // sim.updateTime();


  // Save data in test directory
  strcpy(save.dir, "../TestData/CPUHDF5");
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

  ParallelEnv env(0, NULL, 2, 2, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env, cfl, Ng, gamma);
  SRMHD model(&d);
  Weno3 weno(&d);
  FVS fluxMethod(&d, &weno, &model);
  ParallelOutflow bcs(&d, &env);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveDataHDF5 save(&d, &env, "../TestData/CPUHDF5/RK2SrmhdOutflowOTVSF", ParallelSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // sim.evolve();
  sim.updateTime();
  // sim.updateTime();


  // Save data in test directory
  strcpy(save.dir, "../TestData/CPUHDF5");
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

  ParallelEnv env(0, NULL, 2, 2, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env, cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip);
  SRMHD model(&d);
  Weno3 weno(&d);
  FVS fluxMethod(&d, &weno, &model);
  ParallelOutflow bcs(&d, &env);
  Simulation sim(&d, &env);
  KHRandomInstabilitySingleFluid init(&d, 1, seed);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveDataHDF5 save(&d, &env, "../TestData/CPUHDF5/RK2SrmhdOutflowKHRandomInstabilitySF", ParallelSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/CPUHDF5");
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

  ParallelEnv env(0, NULL, 2, 2, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env, cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip);
  SRMHD model(&d);
  Weno3 weno(&d);
  FVS fluxMethod(&d, &weno, &model);
  ParallelPeriodic bcs(&d, &env);
  Simulation sim(&d, &env);
  KHRandomInstabilitySingleFluid init(&d, 1, seed);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveDataHDF5 save(&d, &env, "../TestData/CPUHDF5/RK2SrmhdPeriodicKHRandomInstabilitySF", ParallelSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/CPUHDF5");
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

  ParallelEnv env(0, NULL, 2, 2, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env, cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip);
  SRMHD model(&d);
  Weno3 weno(&d);
  FVS fluxMethod(&d, &weno, &model);
  ParallelFlow bcs(&d, &env);
  Simulation sim(&d, &env);
  KHRandomInstabilitySingleFluid init(&d, 1, seed);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveDataHDF5 save(&d, &env, "../TestData/CPUHDF5/RK2SrmhdFlowKHRandomInstabilitySF", ParallelSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/CPUHDF5");
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

  ParallelEnv env(0, NULL, 2, 2, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  SRMHD model(&d);
  Weno3 weno(&d);
  FVS fluxMethod(&d, &weno, &model);
  ParallelOutflow bcs(&d, &env);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveDataHDF5 save(&d, &env, "../TestData/CPUHDF5/RK2SrmhdOutflowBrioWuSF", ParallelSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/CPUHDF5");
  strcpy(save.app, "RK2SrmhdOutflowBrioWuSF");

  save.saveAll();
}

TEST(RK2OutputConsistentWithSerial, RK2SrmhdPeriodicBrioWuSF)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  ParallelEnv env(0, NULL, 2, 2, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  SRMHD model(&d);
  Weno3 weno(&d);
  FVS fluxMethod(&d, &weno, &model);
  ParallelPeriodic bcs(&d, &env);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveDataHDF5 save(&d, &env, "../TestData/CPUHDF5/RK2SrmhdPeriodicBrioWuSF", ParallelSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/CPUHDF5");
  strcpy(save.app, "RK2SrmhdPeriodicBrioWuSF");

  save.saveAll();
}

TEST(RK2OutputConsistentWithSerial, RK2SrmhdFlowBrioWuSF)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  ParallelEnv env(0, NULL, 2, 2, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  SRMHD model(&d);
  Weno3 weno(&d);
  FVS fluxMethod(&d, &weno, &model);
  ParallelFlow bcs(&d, &env);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveDataHDF5 save(&d, &env, "../TestData/CPUHDF5/RK2SrmhdFlowBrioWuSF", ParallelSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/CPUHDF5");
  strcpy(save.app, "RK2SrmhdFlowBrioWuSF");

  save.saveAll();
}
#endif



#if 0

// Tests which do not currently pass

TEST(RK2OutputConsistentWithSerial, RK2SrrmhdOutflowOTVortexSingleFluid)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  ParallelEnv env(0, NULL, 2, 2, 1, 1);
  Data d(30, 30, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  SRRMHD model(&d);
  FVS fluxMethod(&d, &model);
  ParallelOutflow bcs(&d, &env);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/CPU");
  strcpy(save.app, "RK2SrrmhdOutflowOTVortexSingleFluid");

  save.saveAll();
}
#endif
