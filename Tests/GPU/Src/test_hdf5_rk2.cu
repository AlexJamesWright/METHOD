#include "gtest/gtest.h"
#include "srrmhd.h"
#include "srmhd.h"
#include "simulation.h"
#include "serialSaveDataHDF5.h"
#include "simData.h"
#include "initFunc.h"
#include "RK2.h"
#include "fluxVectorSplitting.h"
#include <cstdlib>


/*
 Assumptions:
   RKRandomInstabilitySingleFluid is tested in 2D only
   BrioWuSingleFluid is tested in 1D only
*/


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
  SRMHD model(&d);
  FVS fluxMethod(&d, &model);
  Periodic bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveDataHDF5 save(&d, &env, "../TestData/GPUHDF5/RK2SrmhdPeriodicOTVSF", SerialSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // sim.evolve();
  sim.updateTime();
  // sim.updateTime();

  // Save data in test directory
  // This currently needs to be set in the save() function above as well
  strcpy(save.dir, "../TestData/GPUHDF5");
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
  SRMHD model(&d);
  FVS fluxMethod(&d, &model);
  Outflow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveDataHDF5 save(&d, &env, "../TestData/GPUHDF5/RK2SrmhdOutflowOTVSF", SerialSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // sim.evolve();
  sim.updateTime();
  // sim.updateTime();

  // Save data in test directory
  // This currently needs to be set in the save() function above as well
  strcpy(save.dir, "../TestData/GPUHDF5");
  strcpy(save.app, "RK2SrmhdOutflowOTVSF");

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
  SRMHD model(&d);
  FVS fluxMethod(&d, &model);
  Outflow bcs(&d);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);

  SerialSaveDataHDF5 save(&d, &env, "../TestData/GPUHDF5/RK2SrmhdOutflowBrioWuSF", SerialSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();
  //sim.updateTime();
  // sim.updateTime();

  // Save data in test directory
  // This currently needs to be set in the save() function above as well
  strcpy(save.dir, "../TestData/GPUHDF5");
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
  SRMHD model(&d);
  FVS fluxMethod(&d, &model);
  Periodic bcs(&d);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);

  SerialSaveDataHDF5 save(&d, &env, "../TestData/GPUHDF5/RK2SrmhdPeriodicBrioWuSF", SerialSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();
  //sim.updateTime();
  // sim.updateTime();

  // Save data in test directory
  // This currently needs to be set in the save() function above as well
  strcpy(save.dir, "../TestData/GPUHDF5");
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
  SRMHD model(&d);
  FVS fluxMethod(&d, &model);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RK2 timeInt(&d, &model, &bcs, &fluxMethod);

  SerialSaveDataHDF5 save(&d, &env, "../TestData/GPUHDF5/RK2SrmhdFlowBrioWuSF", SerialSaveDataHDF5::OUTPUT_ALL);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  sim.evolve();
  //sim.updateTime();
  // sim.updateTime();

  // Save data in test directory
  // This currently needs to be set in the save() function above as well
  strcpy(save.dir, "../TestData/GPUHDF5");
  strcpy(save.app, "RK2SrmhdFlowBrioWuSF");

  save.saveAll();
}



