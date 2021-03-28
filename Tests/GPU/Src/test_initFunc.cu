#include "gtest/gtest.h"
#include "simData.h"
#include "initFunc.h"
#include "boundaryConds.h"
#include "simulation.h"
#include "srmhd.h"
#include "twoFluidEMHD.h"
#include "serialEnv.h"
#include <iostream>

TEST(InitialFunc, BaseConstructor)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data data(100, 10, 10, 0, 1, -0.5, 0.5, -0.1, 0.1, 0.8, &env);
  SRMHD model(&data);
  Periodic bcs(&data);
  Simulation sim(&data, &env);
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
  SerialEnv env(0, NULL, 1, 1, 1);
  Data data(100, 10, 0, 0, 1, 0, 1, -0.1, 0.1, 0.8, &env);
  SRMHD model(&data);
  Periodic bcs(&data);
  Simulation sim(&data, &env);
  OTVortexSingleFluid init(&data);

  EXPECT_NEAR(data.prims[data.id(0, 0, 0, 0)], 0.2210485321, 0.0000000001);
  EXPECT_NEAR(data.prims[data.id(0, 99, 9, 0)], 0.2210485321, 0.0000000001);
  EXPECT_NEAR(data.prims[data.id(1, 35, 5, 0)], -0.4045084972, 0.0000000001);
  EXPECT_NEAR(data.prims[data.id(2, 34, 2, 0)], 0.4704403845, 0.000000001);
  EXPECT_EQ(data.prims[data.id(3, 50, 5, 0)], 0);
  EXPECT_NEAR(data.prims[data.id(4, 85, 3, 0)], 0.1326291192, 0.0000000001);
  EXPECT_NEAR(data.prims[data.id(5, 33, 12, 0)], 0.2282194806, 0.000000001);
  EXPECT_NEAR(data.prims[data.id(6, 67, 2, 0)], 0.2798703901, 0.0000000001);
  EXPECT_EQ(data.prims[data.id(7, 99, 9, 0)], 0);

}


// Setting up function across different axes should be same after rotation
TEST(InitialFunc, BrioWuTwoFluidFunc)
{
  // Discontinuity in x direction
  SerialEnv env(0, NULL, 1, 1, 1);
  Data dx(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.8, &env);
  TwoFluidEMHD modelx(&dx);
  Periodic bcsx(&dx);
  Simulation simx(&dx, &env);
  BrioWuTwoFluid initx(&dx, 0);



  // Discontinuity in y direction
  SerialEnv env2(0, NULL, 1, 1, 1);
  Data dy(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.8, &env2);
  TwoFluidEMHD modely(&dy);
  Periodic bcsy(&dy);
  Simulation simy(&dy, &env2);
  BrioWuTwoFluid inity(&dy, 1);


  // Discontinuity in z direction
  SerialEnv env3(0, NULL, 1, 1, 1);
  Data dz(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.8, &env3);
  TwoFluidEMHD modelz(&dz);
  Periodic bcsz(&dz);
  Simulation simz(&dz, &env3);
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
