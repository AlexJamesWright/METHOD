#include "gtest/gtest.h"
#include "srrmhd.h"
#include "simulation.h"
#include "simData.h"
#include "initFunc.h"
#include "RK2.h"
#include "fluxVectorSplitting.h"
#include <cstdlib>


TEST(RK2, RK2OutputConsistentWithSerial)
{

  /*
    The following was used to gather data to compare the parallel
     version with. No tests are run in the serial version of this test
  */

  Data d(30, 30, 30, 0, 1, 0, 1, 0, 1, 0.004);
  SRRMHD model(&d);
  FVS fluxMethod(&d, &model);
  Simulation sim(&d);
  OTVortexSingleFluid init(&d);
  Outflow bcs(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  SaveData save(&d);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();


  // Save data in test directory
  strcpy(save.dir, "../TestData/Parallel");
  strcpy(save.app, "RK2");
  save.saveCons();
  save.savePrims();
  save.saveAux();
  save.saveConsts();



}
