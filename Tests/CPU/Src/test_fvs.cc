#include "gtest/gtest.h"
#include "srmhd.h"
#include "twoFluidEMHD.h"
#include "simulation.h"
#include "simData.h"
#include "serialSaveData.h"
#include "initFunc.h"
#include "rkSplit.h"
#include "fluxVectorSplitting.h"
#include "serialEnv.h"
#include <cstdlib>
#include <cmath>
#include <cstdio>

TEST(FVS, SameFnetAsSerial)
/*!
  Determine the flux for the first step of the OTvortex, check it is the same
  as the serial version.
*/
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(20, 20, 20, 0, 1, 0, 1, 0, 1, 0.8, &env);
  SRMHD model(&d);
  Weno3 weno(&d);
  FVS fluxMethod(&d, &weno, &model);
  Outflow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  fluxMethod.F(d.cons, d.prims, d.aux, d.f, d.fnet);

  for (int var(0); var < d.Ncons; var++)
  {
    for (int i(0); i < d.Nx; i++)
    {
      for (int j(0); j < d.Ny; j++)
      {
        for (int k(0); k < d.Nz; k++)
        {
          d.cons[d.id(var, i, j, k)] = d.fnet[d.id(var, i, j, k)];
        }
      }
    }
  }

  // Save data in test directory
  strcpy(save.dir, "../TestData/Serial");
  strcpy(save.app, "FVSFnet");
  save.saveCons();
  save.saveConsts();
}

TEST(FVS, SameXReconstructionAsSerial)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(20, 20, 20, 0, 1, 0, 1, 0, 1, 0.8, &env);
  SRMHD model(&d);
  Weno3 weno(&d);
  FVS fluxMethod(&d, &weno, &model);
  Outflow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  model.fluxVector(d.cons, d.prims, d.aux, d.f, 0);
  fluxMethod.fluxReconstruction(d.cons, d.prims, d.aux, d.f, d.fnet, 0);

  for (int var(0); var < d.Ncons; var++)
  {
    for (int i(0); i < d.Nx-1; i++)
    {
      for (int j(0); j < d.Ny; j++)
      {
        for (int k(0); k < d.Nz; k++)
        {
          d.cons[d.id(var, i, j, k)] = d.fnet[d.id(var, i+1, j, k)]/d.dx - d.fnet[d.id(var, i, j, k)]/d.dx;
        }
      }
    }
  }

  // Save data in test directory
  strcpy(save.dir, "../TestData/Serial");
  strcpy(save.app, "FVSFx");
  save.saveCons();
  save.saveConsts();

}

TEST(FVS, SameYReconstructionAsSerial)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(20, 20, 20, 0, 1, 0, 1, 0, 1, 0.8, &env);
  SRMHD model(&d);
  Weno3 weno(&d);
  FVS fluxMethod(&d, &weno, &model);
  Outflow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  model.fluxVector(d.cons, d.prims, d.aux, d.f, 1);
  fluxMethod.fluxReconstruction(d.cons, d.prims, d.aux, d.f, d.fnet, 1);

  for (int var(0); var < d.Ncons; var++)
  {
    for (int i(0); i < d.Nx; i++)
    {
      for (int j(0); j < d.Ny; j++)
      {
        for (int k(0); k < d.Nz; k++)
        {
          d.cons[d.id(var, i, j, k)] = d.fnet[d.id(var, i, j+1, k)]/d.dy - d.fnet[d.id(var, i, j, k)]/d.dy;
        }
      }
    }
  }

  // Save data in test directory
  strcpy(save.dir, "../TestData/Serial");
  strcpy(save.app, "FVSFy");
  save.saveCons();
  save.saveConsts();

}

TEST(FVS, SameZReconstructionAsSerial)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(20, 20, 20, 0, 1, 0, 1, 0, 1, 0.8, &env);
  SRMHD model(&d);
  Weno3 weno(&d);
  FVS fluxMethod(&d, &weno, &model);
  Outflow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  model.fluxVector(d.cons, d.prims, d.aux, d.f, 2);
  fluxMethod.fluxReconstruction(d.cons, d.prims, d.aux, d.f, d.fnet, 2);

  for (int var(0); var < d.Ncons; var++)
  {
    for (int i(0); i < d.Nx; i++)
    {
      for (int j(0); j < d.Ny; j++)
      {
        for (int k(0); k < d.Nz; k++)
        {
          d.cons[d.id(var, i, j, k)] = d.fnet[d.id(var, i, j, k+1)]/d.dz - d.fnet[d.id(var, i, j, k)]/d.dz;
        }
      }
    }
  }

  // Save data in test directory
  strcpy(save.dir, "../TestData/Serial");
  strcpy(save.app, "FVSFz");
  save.saveCons();
  save.saveConsts();

}
