#include "gtest/gtest.h"
#include "../Src/simulation.h"
#include <cstdlib>

namespace {

  TEST(Simulation, Null)
  {
    Simulation sim(100, 10, 0, 1, -0.5, 0.5, 0.8);
    EXPECT_EQ(sim.Nx, 100);
    EXPECT_EQ(sim.Ny, 10);
    EXPECT_EQ(sim.xmin, 0.0);
    EXPECT_EQ(sim.xmax, 1.0);
    EXPECT_EQ(sim.ymin, -0.5);
    EXPECT_EQ(sim.ymax, 0.5);
    EXPECT_EQ(sim.endTime, 0.8);
    EXPECT_EQ(sim.cfl, 0.5);
    EXPECT_EQ(sim.Ng, 4);
    EXPECT_EQ(sim.gamma, 5.0/3);
    EXPECT_EQ(sim.sigma, 10.0);

  }


}
