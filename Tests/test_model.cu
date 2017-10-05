#include "gtest/gtest.h"
#include "../Src/model.h"
#include "../Src/simulation.h"
#include <cstdlib>

namespace {

  TEST(SRMHD, ModelDefaultConstructor)
  {
    SRMHD model;
    EXPECT_EQ(model.Ncons, 9);
    EXPECT_EQ(model.Nprims, 8);
    EXPECT_EQ(model.Naux, 9);
  }

  TEST(SRMHD, ModelConstructor)
  {
    EXPECT_EQ(1, 1);
  }

}
