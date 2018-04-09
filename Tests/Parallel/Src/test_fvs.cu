#include "gtest/gtest.h"
#include "srmhd.h"
#include "twoFluidEMHD.h"
#include "simulation.h"
#include "simData.h"
#include "initFunc.h"
#include "rkSplit.h"
#include "fluxVectorSplitting.h"
#include <cstdlib>
#include <cmath>
#include <cstdio>

#define ID(variable, idx, jdx, kdx)  ((variable)*(dx.Nx)*(dx.Ny)*(dx.Nz) + (idx)*(dx.Ny)*(dx.Nz) + (jdx)*(dx.Nz) + (kdx))


namespace
{

  TEST(FVS, SameFnetAsSerial)
  {
    Data d(50, 50, 50, 0, 1, 0, 1, 0, 1, 0.8);
    SRMHD model(&d);
    FVS fluxMethod(&d, &model);
    Simulation sim(&d);
    OTVortexSingleFluid init(&d);
    Outflow bcs(&d);
    RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
    SaveData save(&d);
    sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

    fluxMethod.F(d.cons, d.prims, d.aux, d.f, d.fnet);
    EXPECT_NEAR(d.fnet[d.id(0, 19, 19, 19)], -0.0140313284779348, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(1, 19, 19, 19)], 1.3102506615156173, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(2, 19, 19, 19)], -0.5130098084853181, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(3, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(4, 19, 19, 19)], -1.0707540789952823, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(5, 19, 19, 19)], 0.5421499985520821, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(6, 19, 19, 19)], 0.9676495105836604, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(7, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(8, 19, 19, 19)], 0.0000000000000000, 1.0e-16);

  }

  TEST(FVS, SameXReconstructionAsSerial)
  {
    Data d(30, 30, 30, 0, 1, 0, 1, 0, 1, 0.8);
    TwoFluidEMHD model(&d);
    FVS fluxMethod(&d, &model);
    Simulation sim(&d);
    BrioWuTwoFluid init(&d, 0, 0);
    Outflow bcs(&d);
    RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
    SaveData save(&d);
    sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

    fluxMethod.fluxReconstruction(d.cons, d.prims, d.aux, d.f, d.fnet, 0);
    EXPECT_NEAR(d.fnet[d.id(0, 19, 19, 19)], -12.7499999999999325, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(1, 19, 19, 19)], -13.5000000000000036, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(2, 19, 19, 19)], 15.0000000000000071, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(3, 19, 19, 19)], 15.0000000000000071, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(4, 19, 19, 19)], -20.2499999999999716, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(5, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(6, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(7, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(8, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(9, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(10, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(11, 19, 19, 19)], -29.9999999999999964, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(12, 19, 19, 19)], -29.9999999999999964, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(13, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(14, 19, 19, 19)], -30.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(15, 19, 19, 19)], 30.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(16, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(17, 19, 19, 19)], 0.0000000000000000, 1.0e-16);

  }

  TEST(FVS, SameYReconstructionAsSerial)
  {
    Data d(30, 30, 30, 0, 1, 0, 1, 0, 1, 0.8);
    TwoFluidEMHD model(&d);
    FVS fluxMethod(&d, &model);
    Simulation sim(&d);
    BrioWuTwoFluid init(&d, 1, 0);
    Outflow bcs(&d);
    RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
    SaveData save(&d);
    sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

    fluxMethod.fluxReconstruction(d.cons, d.prims, d.aux, d.f, d.fnet, 1);
    EXPECT_NEAR(d.fnet[d.id(0, 19, 19, 19)], -12.7499999999999325, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(1, 19, 19, 19)], 15.0000000000000071, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(2, 19, 19, 19)], -13.5000000000000036, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(3, 19, 19, 19)], 15.0000000000000071, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(4, 19, 19, 19)], -20.2499999999999716, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(5, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(6, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(7, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(8, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(9, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(10, 19, 19, 19)], -29.9999999999999964, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(11, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(12, 19, 19, 19)], -29.9999999999999964, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(13, 19, 19, 19)], 30.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(14, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(15, 19, 19, 19)], -30.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(16, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(17, 19, 19, 19)], 0.0000000000000000, 1.0e-16);

  }

  TEST(FVS, SameZReconstructionAsSerial)
  {
    Data d(30, 30, 30, 0, 1, 0, 1, 0, 1, 0.8);
    TwoFluidEMHD model(&d);
    FVS fluxMethod(&d, &model);
    Simulation sim(&d);
    BrioWuTwoFluid init(&d, 2, 0);
    Outflow bcs(&d);
    RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
    SaveData save(&d);
    sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

    fluxMethod.fluxReconstruction(d.cons, d.prims, d.aux, d.f, d.fnet, 2);
    EXPECT_NEAR(d.fnet[d.id(0, 19, 19, 19)], -12.7499999999999325, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(1, 19, 19, 19)], 15.0000000000000071, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(2, 19, 19, 19)], 15.0000000000000071, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(3, 19, 19, 19)], -13.5000000000000036, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(4, 19, 19, 19)], -20.2499999999999716, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(5, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(6, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(7, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(8, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(9, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(10, 19, 19, 19)], -29.9999999999999964, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(11, 19, 19, 19)], -29.9999999999999964, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(12, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(13, 19, 19, 19)], -30.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(14, 19, 19, 19)], 30.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(15, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(16, 19, 19, 19)], 0.0000000000000000, 1.0e-16);
    EXPECT_NEAR(d.fnet[d.id(17, 19, 19, 19)], 0.0000000000000000, 1.0e-16);

  }

}
