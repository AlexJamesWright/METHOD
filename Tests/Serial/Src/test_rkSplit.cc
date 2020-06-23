#include "gtest/gtest.h"
#include "srrmhd.h"
#include "srmhd.h"
#include "simulation.h"
#include "simData.h"
#include "saveData.h"
#include "serialSaveData.h"
#include "initFunc.h"
#include "rkSplit.h"
#include "fluxVectorSplitting.h"
#include "platformEnv.h"
#include <cstdlib>

#if 1
// RKRandomInstabilitySingleFluid
TEST(RK2OutputConsistentWithSerial, RKSplitSrmhdOutflowKHRandomInstabilitySF)
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

  PlatformEnv env(0, NULL, 1, 1, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env, cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip);
  SRMHD model(&d);
  FVS fluxMethod(&d, &model);
  Outflow bcs(&d);
  Simulation sim(&d, &env);
  KHRandomInstabilitySingleFluid init(&d, 1, seed);
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/Serial");
  strcpy(save.app, "RKSplitSrmhdOutflowKHRandomInstabilitySF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}

TEST(RK2OutputConsistentWithSerial, RKSplitSrmhdPeriodicKHRandomInstabilitySF)
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

  PlatformEnv env(0, NULL, 1, 1, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env, cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip);
  SRMHD model(&d);
  FVS fluxMethod(&d, &model);
  Periodic bcs(&d);
  Simulation sim(&d, &env);
  KHRandomInstabilitySingleFluid init(&d, 1, seed);
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/Serial");
  strcpy(save.app, "RKSplitSrmhdPeriodicKHRandomInstabilitySF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}

TEST(RK2OutputConsistentWithSerial, RKSplitSrmhdFlowKHRandomInstabilitySF)
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

  PlatformEnv env(0, NULL, 1, 1, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env, cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip);
  SRMHD model(&d);
  FVS fluxMethod(&d, &model);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  KHRandomInstabilitySingleFluid init(&d, 1, seed);
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/Serial");
  strcpy(save.app, "RKSplitSrmhdFlowKHRandomInstabilitySF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}

// BrioWuSingleFluid

TEST(RKSplitOutputConsistentWithSerial, RKSplitSrmhdOutflowBrioWuSF)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  PlatformEnv env(0, NULL, 1, 1, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  SRMHD model(&d);
  FVS fluxMethod(&d, &model);
  Outflow bcs(&d);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/Serial");
  strcpy(save.app, "RKSplitSrmhdOutflowBrioWuSF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}

TEST(RKSplitOutputConsistentWithSerial, RKSplitSrmhdPeriodicBrioWuSF)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  PlatformEnv env(0, NULL, 1, 1, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  SRMHD model(&d);
  FVS fluxMethod(&d, &model);
  Periodic bcs(&d);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/Serial");
  strcpy(save.app, "RKSplitSrmhdPeriodicBrioWuSF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}

TEST(RKSplitOutputConsistentWithSerial, RKSplitSrmhdFlowBrioWuSF)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  PlatformEnv env(0, NULL, 1, 1, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  SRMHD model(&d);
  FVS fluxMethod(&d, &model);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/Serial");
  strcpy(save.app, "RKSplitSrmhdFlowBrioWuSF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}
#endif

#if 0
TEST(RKSplitOutputConsistentWithSerial, RKSplitSrrmhdOutflowOTVortexSingleFluidRKSplit)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  PlatformEnv env(0, NULL, 1, 1, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  SRRMHD model(&d);
  FVS fluxMethod(&d, &model);
  Outflow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/Serial");
  strcpy(save.app, "RKSplitSrrmhdOutflowOTVortexSingleFluid");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}
#endif 
