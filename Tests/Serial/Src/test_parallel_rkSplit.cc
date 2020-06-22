#include "gtest/gtest.h"
#include "srrmhd.h"
#include "srmhd.h"
#include "simulation.h"
#include "simData.h"
#include "parallelSaveData.h"
#include "initFunc.h"
#include "rkSplit.h"
#include "fluxVectorSplitting.h"
#include "platformEnv.h"
#include <cstdlib>

#if 1
TEST(RKSplitOutputConsistentWithSerial, RKSplitSrmhdOutflowBrioWuSF)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  PlatformEnv env(0, NULL, 3, 1, 1, 1);
  Data d(30, 0, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  SRMHD model(&d);
  FVS fluxMethod(&d, &model);
  ParallelOutflow bcs(&d, &env);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/ParallelMPI");
  strcpy(save.app, "RKSplitSrmhdOutflowBrioWuSF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}
#endif

#if 0
TEST(RKSplitOutputConsistentWithSerial, RKSplitSrrmhdOutflowOTVortexSingleFluid)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  PlatformEnv env(0, NULL, 2, 2, 1, 1);
  Data d(30, 30, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  SRRMHD model(&d);
  FVS fluxMethod(&d, &model);
  ParallelOutflow bcs(&d, &env);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  ParallelSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/ParallelMPI");
  strcpy(save.app, "RKSplitSrrmhdOutflowOTVortexSingleFluid");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}
#endif 
