#include "gtest/gtest.h"
#include "simulation.h"
#include "simData.h"
#include "cudaErrorCheck.h"
#include <cstdlib>

namespace {

  Data data(100, 10, 0, 1, -0.5, 0.5, 0.8);

  TEST(Simulation, dataInitialisation)
  {
    Simulation sim(&data);

    // Check standard data
    EXPECT_EQ(sim.data->Nx, 100);
    EXPECT_EQ(sim.data->Ny, 10);
    EXPECT_EQ(sim.data->xmin, 0.0);
    EXPECT_EQ(sim.data->xmax, 1.0);
    EXPECT_EQ(sim.data->ymin, -0.5);
    EXPECT_EQ(sim.data->ymax, 0.5);
    EXPECT_EQ(sim.data->endTime, 0.8);
    EXPECT_EQ(sim.data->cfl, 0.5);
    EXPECT_EQ(sim.data->Ng, 4);
    EXPECT_EQ(sim.data->gamma, 5.0/3);
    EXPECT_EQ(sim.data->sigma, 0.0);

    // Check domain
    /* Note: Ng = 4 */
    EXPECT_EQ(sim.data->x[4], 0.005);
    EXPECT_EQ(sim.data->x[4+1], 0.015);
    EXPECT_EQ(sim.data->y[4], -0.45);
    EXPECT_EQ(sim.data->y[4+1], -0.35);

  }


}
