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

      Data d(30, 30, 30, 0, 1, 0, 1, 0, 1, 0.8);
      TwoFluidEMHD model(&d);
      FVS fluxMethod(&d, &model);
      Simulation sim(&d);
      BrioWuTwoFluid init(&d, 0, 0);
      Outflow bcs(&d);
      RK2 timeInt(&d, &model, &bcs, &fluxMethod);
      SaveData save(&d);
      sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

      sim.updateTime();

      printf("Results:\n");
      for (int var(0); var < d.Ncons; var++) {
        printf("%19.16f\n", d.cons[d.id(var, 19, 0, 19)]);
      }

    }
}
