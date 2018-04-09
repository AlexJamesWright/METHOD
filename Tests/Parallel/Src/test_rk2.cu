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

      EXPECT_NEAR(d.cons[d.id(0, 19, 0, 19)], 0.1619806802286889, 1.0e-16);
      EXPECT_NEAR(d.cons[d.id(1, 19, 0, 19)], 0.0125899405040632, 1.0e-16);
      EXPECT_NEAR(d.cons[d.id(2, 19, 0, 19)], -0.0143295913904252, 1.0e-16);
      EXPECT_NEAR(d.cons[d.id(3, 19, 0, 19)], -0.0143295913904252, 1.0e-16);
      EXPECT_NEAR(d.cons[d.id(4, 19, 0, 19)], 1.2942058736254853, 1.0e-16);
      EXPECT_NEAR(d.cons[d.id(5, 19, 0, 19)], 0.0000000000000000, 1.0e-16);
      EXPECT_NEAR(d.cons[d.id(6, 19, 0, 19)], 0.0000000000000000, 1.0e-16);
      EXPECT_NEAR(d.cons[d.id(7, 19, 0, 19)], 0.0000000000000000, 1.0e-16);
      EXPECT_NEAR(d.cons[d.id(8, 19, 0, 19)], 0.0000000000000000, 1.0e-16);
      EXPECT_NEAR(d.cons[d.id(9, 19, 0, 19)], 0.0000000000000000, 1.0e-16);
      EXPECT_NEAR(d.cons[d.id(10, 19, 0, 19)], 0.5000000000000000, 1.0e-16);
      EXPECT_NEAR(d.cons[d.id(11, 19, 0, 19)], -0.9713408172191497, 1.0e-16);
      EXPECT_NEAR(d.cons[d.id(12, 19, 0, 19)], -0.9713408172191497, 1.0e-16);
      EXPECT_NEAR(d.cons[d.id(13, 19, 0, 19)], 0.0000000000000000, 1.0e-16);
      EXPECT_NEAR(d.cons[d.id(14, 19, 0, 19)], 0.0286591827808504, 1.0e-16);
      EXPECT_NEAR(d.cons[d.id(15, 19, 0, 19)], -0.0286591827808504, 1.0e-16);
      EXPECT_NEAR(d.cons[d.id(16, 19, 0, 19)], 0.0000000000000000, 1.0e-16);
      EXPECT_NEAR(d.cons[d.id(17, 19, 0, 19)], 0.0000000000000000, 1.0e-16);

    }
}
