#include "gtest/gtest.h"
#include "srrmhd.h"
#include "simulation.h"
#include "simData.h"
#include "boundaryConds.h"
#include "initFunc.h"
#include "fluxVectorSplitting.h"
#include "serialEnv.h"
#include <cstdlib>
#include <cmath>
#include <stdio.h>


/* ######################### Test model constructor ########################*/

TEST(SRRMHD, Constructor)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(100, 10, 0, 0, 1, -0.5, 0.5, -0.1, 0.1, 0.8, &env);
  SRRMHD model(&d);
  EXPECT_EQ(d.Ncons, 14);
  EXPECT_EQ(d.Nprims, 11);
  EXPECT_EQ(d.Naux, 17);

}



/* ######################### Test flux vector splitting ########################*/

TEST(SRRMHD, FluxVectorSplittingStationary)
{
  double tol(1.0e-15);
  // Set up
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 1.0, &env, 0.5, 4, 5.0/3.0, 1000.0, 0.5);
  SRRMHD model(&d);
  FVS fluxMethod(&d, &model);
  Periodic bcs(&d);
  Simulation sim(&d, &env);

  // Set state to stationary equilibrium state
  for (int i(0); i < d.Nx; i++) {
    for (int j(0); j < d.Ny; j++) {
      for (int k(0); k < d.Nz; k++) {
        d.prims[d.id(0, i, j, k)] = 0.5; // Require non-zero density
        d.prims[d.id(1, i, j, k)] = 0.0;
        d.prims[d.id(2, i, j, k)] = 0.0;
        d.prims[d.id(3, i, j, k)] = 0.0;
        d.prims[d.id(4, i, j, k)] = 0.0;
        d.prims[d.id(5, i, j, k)] = 0.0;
        d.prims[d.id(6, i, j, k)] = 0.0;
        d.prims[d.id(7, i, j, k)] = 0.0;
        d.prims[d.id(8, i, j, k)] = 0.0;
        d.prims[d.id(9, i, j, k)] = 0.0;
        d.prims[d.id(10, i, j, k)] = 0.0;
      }
    }
  }

  model.primsToAll(d.cons, d.prims, d.aux);

  // System is stationary, there should be zero flux
  // x-direction
  fluxMethod.fluxReconstruction(d.cons, d.prims, d.aux, d.f, d.fnet, 0);
  for (int i(d.Ng); i < d.Nx-d.Ng; i++) {
    for (int j(d.Ng); j < d.Ny-d.Ng; j++) {
      for (int k(d.Ng); k < d.Nz-d.Ng; k++) {
        for (int var(0); var < d.Ncons; var++) {
          EXPECT_NEAR(d.fnet[d.id(var, i, j, k)], 0.0, tol);
        }
      }
    }
  }
  // y-direction
  fluxMethod.fluxReconstruction(d.cons, d.prims, d.aux, d.f, d.fnet, 0);
  for (int i(d.Ng); i < d.Nx-d.Ng; i++) {
    for (int j(d.Ng); j < d.Ny-d.Ng; j++) {
      for (int k(d.Ng); k < d.Nz-d.Ng; k++) {
        for (int var(0); var < d.Ncons; var++) {
          EXPECT_NEAR(d.fnet[d.id(var, i, j, k)], 0.0, tol);
        }
      }
    }
  }
  // z-direction
  fluxMethod.fluxReconstruction(d.cons, d.prims, d.aux, d.f, d.fnet, 0);
  for (int i(d.Ng); i < d.Nx-d.Ng; i++) {
    for (int j(d.Ng); j < d.Ny-d.Ng; j++) {
      for (int k(d.Ng); k < d.Nz-d.Ng; k++) {
        for (int var(0); var < d.Ncons; var++) {
          EXPECT_NEAR(d.fnet[d.id(var, i, j, k)], 0.0, tol);
        }
      }
    }
  }
}


/* ######################### Test getPrimVars transform ########################*/

TEST(SRRMHD, Prims2Cons2Prims)
{
  const double tol = 1.49011612e-8;   // Tolerance of rootfinder
  SerialEnv env(0, NULL, 1, 1, 1);
  SerialEnv env2(0, NULL, 1, 1, 1);
  Data d(10, 10, 0, 0, 1, 0, 1, 0, 1, 1.0, &env);
  Data d2(10, 10, 0, 0, 1, 0, 1, 0, 1, 1.0, &env2);
  SRRMHD model(&d);
  SRRMHD model2(&d2);
  Periodic bcs(&d);
  Periodic bcs2(&d2);
  Simulation sim(&d, &env);
  Simulation sim2(&d2, &env2);
  OTVortexSingleFluid init(&d);
  OTVortexSingleFluid init2(&d2);

  model2.primsToAll(d2.cons, d2.prims, d2.aux);
  model.primsToAll(d.cons, d.prims, d.aux);

  model2.getPrimitiveVars(d2.cons, d2.prims, d2.aux);


  for (int var(0); var < d.Nprims; var++) {
    for (int i(2); i < d.Nx-2; i++) {
      for (int j(2); j < d.Ny-2; j++) {
        EXPECT_NEAR(d.prims[d.id(var, i, j, 0)], d2.prims[d.id(var, i, j, 0)], tol);
      }
    }
  }
  for (int var(0); var < d.Naux; var++) {
    for (int i(2); i < d.Nx-2; i++) {
      for (int j(2); j < d.Ny-2; j++) {
        EXPECT_NEAR(d.aux[d.id(var, i, j, 0)], d2.aux[d.id(var, i, j, 0)], tol);
      }
    }
  }

  // Set all d2 prims slightly off so rootfind has to do something
  for (int var(0); var < d.Nprims; var++) {
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        d2.prims[d.id(var, i, j, 0)] *= 0.9;
      }
    }
  }
  for (int var(0); var < d.Naux; var++) {
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        d2.aux[d.id(var, i, j, 0)] *= 0.9;
      }
    }
  }
  // Solve and re-check
  model2.getPrimitiveVars(d2.cons, d2.prims, d2.aux);

  for (int var(0); var < d.Nprims; var++) {
    for (int i(2); i < d.Nx-2; i++) {
      for (int j(2); j < d.Ny-2; j++) {
        EXPECT_NEAR(d.prims[d.id(var, i, j, 0)], d2.prims[d.id(var, i, j, 0)], tol);
      }
    }
  }
  for (int var(0); var < d.Naux; var++) {
    for (int i(2); i < d.Nx-2; i++) {
      for (int j(2); j < d.Ny-2; j++) {
        EXPECT_NEAR(d.aux[d.id(var, i, j, 0)], d2.aux[d.id(var, i, j, 0)], tol);
      }
    }
  }
}
