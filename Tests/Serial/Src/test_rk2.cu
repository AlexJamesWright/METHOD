#include "gtest/gtest.h"
#include "twoFluidEMHD.h"
#include "simulation.h"
#include "simData.h"
#include "initFunc.h"
#include "RK2.h"
#include "fluxVectorSplitting.h"
#include <cstdlib>
#include <cmath>
#include <cstdio>

namespace
{

  TEST(RK2, RK2OutputConsistentWithSerial)
    {

      /*
        The following was used to gather data to compare the parallel
         version with. No tests are run in the serial version of this test
      */

      // Data d(30, 30, 30, 0, 1, 0, 1, 0, 1, 0.8);
      // TwoFluidEMHD model(&d);
      // FVS fluxMethod(&d, &model);
      // Simulation sim(&d);
      // BrioWuTwoFluid init(&d, 0, 0);
      // Outflow bcs(&d);
      // RK2 timeInt(&d, &model, &bcs, &fluxMethod);
      // SaveData save(&d);
      // sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
      //
      // sim.updateTime();



    }
}
