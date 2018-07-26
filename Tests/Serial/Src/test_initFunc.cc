#include "gtest/gtest.h"
#include "simData.h"
#include "initFunc.h"
#include "simulation.h"
#include "srmhd.h"
#include "twoFluidEMHD.h"
#include <iostream>

namespace
{

  TEST(InitialFunc, baseConstructor)
  {
    Data data(100, 10, 10, 0, 1, -0.5, 0.5, -0.1, 0.1, 0.8);
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

  TEST(InitialFunc, OTVortexSingleFluidFunc)
  {
    Data data(100, 10, 2, 0, 1, 0, 1, -0.1, 0.1, 0.8);
    SRMHD model(&data);
    Simulation sim(&data);
    OTVortexSingleFluid init(&data);

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


  // Setting up function across different axes should be same after rotation
  TEST(InitialFunc, BrioWuTwoFluidFunc)
  {
    // Discontinuity in x direction
    Data dx(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.8);
    TwoFluidEMHD modelx(&dx);
    Simulation simx(&dx);
    BrioWuTwoFluid initx(&dx, 0);



    // Discontinuity in y direction
    Data dy(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.8);
    TwoFluidEMHD modely(&dy);
    Simulation simy(&dy);
    BrioWuTwoFluid inity(&dy, 1);


    // Discontinuity in z direction
    Data dz(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.8);
    TwoFluidEMHD modelz(&dz);
    Simulation simz(&dz);
    BrioWuTwoFluid initz(&dz, 2);

    for (int var(0); var < dx.Ncons; var++) {
      for (int i(0); i < dx.Nx; i++) {
        for (int j(0); j < dy.Ny; j++) {
          for (int k(0); k < dz.Nz; k++) {
            EXPECT_EQ(dx.cons[dx.id(var, i, j, k)], dy.cons[dy.id(var, j, i, k)]);
            EXPECT_EQ(dx.cons[dx.id(var, i, j, k)], dz.cons[dy.id(var, k, j, i)]);
            EXPECT_EQ(dy.cons[dx.id(var, i, j, k)], dz.cons[dy.id(var, i, k, j)]);
          }
        }
      }
    }
  }


}
