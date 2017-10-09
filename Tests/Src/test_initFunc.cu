#include "gtest/gtest.h"
#include "simData.h"
#include "initFunc.h"
#include "simulation.h"

namespace {

  Data data(100, 10, 0, 1, -0.5, 0.5, 0.8);
  SRMHD model(&data);

  Simulation sim(&data);


  TEST(InitialFunc, baseConstructor) {

    InitialFunc init(&data);

    EXPECT_EQ(data.prims[0], 0);
    EXPECT_EQ(data.prims[5], 0);
    EXPECT_EQ(data.prims[data.id(0, 90, 9)], 0);
    EXPECT_EQ(data.prims[data.id(0, 99, 3)], 0);
    EXPECT_EQ(data.prims[data.id(7, 36, 3)], 0);

    EXPECT_EQ(data.cons[0], 0);
    EXPECT_EQ(data.cons[5], 0);
    EXPECT_EQ(data.cons[data.id(3, 90, 9)], 0);
    EXPECT_EQ(data.cons[data.id(2, 99, 3)], 0);
    EXPECT_EQ(data.cons[data.id(8, 36, 3)], 0);

    EXPECT_EQ(data.f[3], 0);
    EXPECT_EQ(data.f[46], 0);
    EXPECT_EQ(data.f[data.id(1, 90, 9)], 0);
    EXPECT_EQ(data.f[data.id(0, 9, 8)], 0);
    EXPECT_EQ(data.f[data.id(6, 77, 5)], 0);

    EXPECT_EQ(data.fnet[0], 0);
    EXPECT_EQ(data.fnet[5], 0);
    EXPECT_EQ(data.fnet[data.id(0, 90, 9)], 0);
    EXPECT_EQ(data.fnet[data.id(0, 99, 3)], 0);
    EXPECT_EQ(data.fnet[data.id(7, 36, 3)], 0);
  }


}
