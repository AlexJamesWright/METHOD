#include "gtest/gtest.h"
#include "srrmhd.h"
#include "srmhd.h"
#include "simulation.h"
#include "simData.h"
#include "parallelSaveData.h"
#include "initFunc.h"
#include "RK2.h"
#include "fluxVectorSplitting.h"
#include "platformEnv.h"
#include <cstdlib>

#if 1
TEST(RK2OutputConsistentWithSerial, RK2SrmhdOutflowBrioWuSF)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  PlatformEnv env(0, NULL, 4, 1, 1, 1);
  Data d(40, 0, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
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
  strcpy(save.dir, "../TestData/ParallelMPI");
  strcpy(save.app, "RK2SrmhdOutflowBrioWuSF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}
#endif

#if 0
TEST(RK2OutputConsistentWithSerial, RK2SrrmhdOutflowOTVortexSingleFluid)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  PlatformEnv env(0, NULL, 2, 2, 1, 1);
  Data d(40, 40, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
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
  strcpy(save.dir, "../TestData/ParallelMPI");
  strcpy(save.app, "RK2SrrmhdOutflowOTVortexSingleFluid");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}
#endif 
