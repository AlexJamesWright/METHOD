#include "gtest/gtest.h"
#include "srrmhd.h"
#include "srmhd.h"
#include "simulation.h"
#include "simData.h"
#include "saveData.h"
#include "serialSaveData.h"
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

  PlatformEnv env(0, NULL, 1, 1, 1, 1);
  Data d(30, 0, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  SRMHD model(&d);
  FVS fluxMethod(&d, &model);
  Outflow bcs(&d);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/Serial");
  strcpy(save.app, "RK2SrmhdOutflowBrioWuSF");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}
#endif

#if 0
TEST(RK2OutputConsistentWithSerial, RK2SrrmhdOutflowOTVortexSingleFluidRK2)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  PlatformEnv env(0, NULL, 1, 1, 1, 1);
  Data d(30, 30, 0, 0, 1, 0, 1, 0, 1, 0.004, &env);
  SRRMHD model(&d);
  FVS fluxMethod(&d, &model);
  Outflow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/Serial");
  strcpy(save.app, "RK2SrrmhdOutflowOTVortexSingleFluid");

  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();
}
#endif 
