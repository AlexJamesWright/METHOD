#include "gtest/gtest.h"
#include "srrmhd.h"
#include "srmhd.h"
#include "simulation.h"
#include "simData.h"
#include "parallelSaveData.h"
#include "parallelBoundaryConds.h"
#include "initFunc.h"
#include "rkSplit.h"
#include "fluxVectorSplitting.h"
#include "parallelEnv.h"
#include <cstdlib>

/*
 Assumptions:
   RKRandomInstabilitySingleFluid is tested in 2D only
   BrioWuSingleFluid is tested in 1D only
*/


// RKOTVSingleFluid
TEST(RKSplitOutputConsistentWithSerial, RKSplitSrmhdPeriodicOTVSF)
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
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  ParallelPeriodic bcs(&d, &env);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/CPU");
  strcpy(save.app, "RKSplitSrmhdPeriodicOTVSF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}

// RKOTVSingleFluid
TEST(RKSplitOutputConsistentWithSerial, RKSplitSrmhdOutflowOTVSF)
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
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  ParallelOutflow bcs(&d, &env);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/CPU");
  strcpy(save.app, "RKSplitSrmhdOutflowOTVSF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}

// RKRandomInstabilitySingleFluid
TEST(RKSplitOutputConsistentWithSerial, RKSplitSrmhdOutflowKHRandomInstabilitySF)
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
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/CPU");
  strcpy(save.app, "RKSplitSrmhdOutflowKHRandomInstabilitySF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}

TEST(RKSplitOutputConsistentWithSerial, RKSplitSrmhdPeriodicKHRandomInstabilitySF)
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
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/CPU");
  strcpy(save.app, "RKSplitSrmhdPeriodicKHRandomInstabilitySF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}

TEST(RKSplitOutputConsistentWithSerial, RKSplitSrmhdFlowKHRandomInstabilitySF)
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
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/CPU");
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

  ParallelEnv env(0, NULL, 2, 2, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  SRMHD model(&d);
  Weno3 weno(&d);
  FVS fluxMethod(&d, &weno, &model);
  ParallelOutflow bcs(&d, &env);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/CPU");
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

  ParallelEnv env(0, NULL, 2, 2, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  SRMHD model(&d);
  Weno3 weno(&d);
  FVS fluxMethod(&d, &weno, &model);
  ParallelPeriodic bcs(&d, &env);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/CPU");
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

  ParallelEnv env(0, NULL, 2, 2, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  SRMHD model(&d);
  Weno3 weno(&d);
  FVS fluxMethod(&d, &weno, &model);
  ParallelFlow bcs(&d, &env);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/CPU");
  strcpy(save.app, "RKSplitSrmhdFlowBrioWuSF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}
