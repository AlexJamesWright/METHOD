#include "gtest/gtest.h"
#include "srrmhd.h"
#include "srmhd.h"
#include "simulation.h"
#include "simData.h"
#include "parallelSaveData.h"
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
  FVS fluxMethod(&d, &model);
  ParallelPeriodic bcs(&d, &env);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // sim.evolve();
  sim.updateTime();
  // sim.updateTime();


  // Save data in test directory
  strcpy(save.dir, "../TestData/MPIGPU");
  strcpy(save.app, "RK2SrmhdPeriodicOTVSF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
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
  FVS fluxMethod(&d, &model);
  ParallelOutflow bcs(&d, &env);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // sim.evolve();
  sim.updateTime();
  // sim.updateTime();


  // Save data in test directory
  strcpy(save.dir, "../TestData/MPIGPU");
  strcpy(save.app, "RK2SrmhdOutflowOTVSF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
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
  FVS fluxMethod(&d, &model);
  ParallelOutflow bcs(&d, &env);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/MPIGPU");
  strcpy(save.app, "RK2SrmhdOutflowBrioWuSF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
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
  FVS fluxMethod(&d, &model);
  ParallelPeriodic bcs(&d, &env);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/MPIGPU");
  strcpy(save.app, "RK2SrmhdPeriodicBrioWuSF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
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
  FVS fluxMethod(&d, &model);
  ParallelFlow bcs(&d, &env);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/MPIGPU");
  strcpy(save.app, "RK2SrmhdFlowBrioWuSF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
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
  strcpy(save.dir, "../TestData/MPIGPU");
  strcpy(save.app, "RK2SrrmhdOutflowOTVortexSingleFluid");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}
#endif
