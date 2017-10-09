#include "gtest/gtest.h"
#include "model.h"
#include "simulation.h"
#include "simData.h"
#include <cstdlib>

namespace {

  TEST(SRMHD, ModelDefaultConstructor)
  {
    SRMHD model;
    EXPECT_EQ(model.Ncons, 9);
    EXPECT_EQ(model.Nprims, 8);
    EXPECT_EQ(model.Naux, 10);
  }


  TEST(SRMHD, ModelConstructor)
  {
    Data data(100, 10, 0, 1, -0.5, 0.5, 0.8);
    SRMHD model(&data);
    EXPECT_EQ(model.data->Ncons, 9);
    EXPECT_EQ(model.data->Nprims, 8);
    EXPECT_EQ(model.data->Naux, 10);

  }

}
