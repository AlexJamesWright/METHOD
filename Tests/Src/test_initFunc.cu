#include "gtest/gtest.h"
#include "simData.h"
#include "initFunc.h"
#include "simulation.h"
#include "srmhd.h"

namespace
{

  TEST(InitialFunc, baseConstructor)
  {
    Data data(100, 10, 1, 0, 1, -0.5, 0.5, -0.1, 0.1, 0.8);
    SRMHD model(&data);
    Simulation sim(&data);
    InitialFunc init(&data);

    EXPECT_EQ(data.prims[0], 0);
    EXPECT_EQ(data.prims[5], 0);
    EXPECT_EQ(data.prims[data.id(0, 90, 9, 7)], 0);
    EXPECT_EQ(data.prims[data.id(0, 99, 3, 7)], 0);
    EXPECT_EQ(data.prims[data.id(7, 36, 3, 7)], 0);

    EXPECT_EQ(data.cons[0], 0);
    EXPECT_EQ(data.cons[5], 0);
    EXPECT_EQ(data.cons[data.id(3, 90, 9, 7)], 0);
    EXPECT_EQ(data.cons[data.id(2, 99, 3, 7)], 0);
    EXPECT_EQ(data.cons[data.id(8, 36, 3, 7)], 0);

    EXPECT_EQ(data.f[3], 0);
    EXPECT_EQ(data.f[46], 0);
    EXPECT_EQ(data.f[data.id(1, 90, 9, 7)], 0);
    EXPECT_EQ(data.f[data.id(0, 9, 8, 7)], 0);
    EXPECT_EQ(data.f[data.id(6, 77, 5, 7)], 0);

    EXPECT_EQ(data.fnet[0], 0);
    EXPECT_EQ(data.fnet[5], 0);
    EXPECT_EQ(data.fnet[data.id(0, 90, 9, 7)], 0);
    EXPECT_EQ(data.fnet[data.id(0, 99, 3, 7)], 0);
    EXPECT_EQ(data.fnet[data.id(7, 36, 3, 7)], 0);
  }

  TEST(InitialFunc, userDefinedConstructor)
  {
    Data data(100, 10, 1, 0, 1, 0, 1, -0.1, 0.1, 0.8);
    SRMHD model(&data);
    Simulation sim(&data);
    OTVortex init(&data);

    EXPECT_NEAR(data.prims[data.id(0, 0, 0, 7)], 0.2210485321, 0.0000000001);
    EXPECT_NEAR(data.prims[data.id(0, 99, 9, 7)], 0.2210485321, 0.0000000001);
    EXPECT_NEAR(data.prims[data.id(1, 35, 5, 7)], -0.4045084972, 0.0000000001);
    EXPECT_NEAR(data.prims[data.id(2, 34, 2, 7)], 0.4704403845, 0.000000001);
    EXPECT_EQ(data.prims[data.id(3, 50, 5, 7)], 0);
    EXPECT_NEAR(data.prims[data.id(4, 85, 3, 7)], 0.1326291192, 0.0000000001);
    EXPECT_NEAR(data.prims[data.id(5, 33, 12, 7)], 0.2282194806, 0.000000001);
    EXPECT_NEAR(data.prims[data.id(6, 67, 2, 7)], 0.2798703901, 0.0000000001);
    EXPECT_EQ(data.prims[data.id(7, 99, 9, 7)], 0);

  }

}
