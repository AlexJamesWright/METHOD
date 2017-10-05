#include "gtest/gtest.h"
#include "../Src/model.h"
#include "../Src/simulation.h"
#include <cstdlib>
#include <iostream>

namespace {

  Simulation sim(100, 100, 0, 1, 0, 1, 1);
  std::cout << "sim.Nx = <<";

  TEST(SRMHD, ModelConstructor)
  {
    SRMHD model(&sim);
    EXPECT_EQ(sim.Ncons, 9);

  }


}
