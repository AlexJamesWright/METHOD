#include "gtest/gtest.h"
#include "srmhd.h"
#include "boundaryConds.h"
#include "simulation.h"
#include "simData.h"
#include "initFunc.h"
#include "fluxVectorSplitting.h"
#include "serialEnv.h"
#include <cstdlib>
#include <cmath>
#include <stdio.h>

/* ######################### Test model constructor ########################*/

TEST(SRMHD, Constructor)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(100, 10, 0, 0, 1, -0.5, 0.5, -0.1, 0.1, 0.8, &env);
  SRMHD model(&d);
  EXPECT_EQ(d.Ncons, 9);
  EXPECT_EQ(d.Nprims, 8);
  EXPECT_EQ(d.Naux, 13);

}



/* ######################### Test flux vector splitting ########################*/

TEST(SRMHD, FluxVectorSplittingStationary)
{
  double tol(1.0e-15);
  // Set up
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 1.0, &env, 0.5, 4, 5.0/3.0, 1000.0, 0.5);
  SRMHD model(&d);
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


/* ######################### Test source contribution ########################*/

TEST(SRMHD, SourceTerm)
{

  // Set up
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 1.0, &env, 0.5, 4, 5.0/3.0, 1000.0, 0.5);
  SRMHD model(&d);
  Periodic bcs(&d);
  Simulation sim(&d, &env);

  // Set cons to something
  for (int i(0); i < d.Nx; i++) {
    for (int j(0); j < d.Ny; j++) {
      for (int k(0); k < d.Nz; k++) {
        for (int var(0); var < d.Ncons; var++) {
          d.cons[d.id(var, i, j, k)] = 3.1415926;
        }
      }
    }
  }
  // Determine source
  model.sourceTerm(d.cons, d.prims, d.aux, d.source);
  for (int i(0); i < d.Nx; i++) {
    for (int j(0); j < d.Ny; j++) {
      for (int k(0); k < d.Nz; k++) {
        for (int var(0); var < d.Ncons; var++) {
          if (var == 8) EXPECT_EQ(d.source[d.id(var, i, j, k)], -3.1415926 / (0.25));
          else EXPECT_EQ(d.source[d.id(var, i, j, k)], 0);
        }
      }
    }
  }
}



/* ######################### Test getPrimVars transform ########################*/

TEST(SRMHD, Prims2Cons2Prims)
{
  const double tol = 1.49011612e-8;   // Tolerance of rootfinder
  SerialEnv env(0, NULL, 1, 1, 1);
  SerialEnv env2(0, NULL, 1, 1, 1);
  Data d(10, 10, 0, 0, 1, 0, 1, 0, 1, 1.0, &env);
  Data d2(10, 10, 0, 0, 1, 0, 1, 0, 1, 1.0, &env2);
  SRMHD model(&d);
  SRMHD model2(&d2);
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
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        for (int k(0); k < d.Nz; k++) {
          EXPECT_NEAR(d.prims[d.id(var, i, j, k)], d2.prims[d.id(var, i, j, k)], tol);
        }
      }
    }
  }
  for (int var(0); var < d.Naux; var++) {
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        for (int k(0); k < d.Nz; k++) {
          EXPECT_NEAR(d.aux[d.id(var, i, j, k)], d2.aux[d.id(var, i, j, k)], tol);
        }
      }
    }
  }

  // Set all d2 prims slightly off so rootfind has to do something
  for (int var(0); var < d.Nprims; var++) {
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        for (int k(0); k < d.Nz; k++) {
          d2.prims[d.id(var, i, j, k)] *= 0.9;
        }
      }
    }
  }
  for (int var(0); var < d.Naux; var++) {
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        for (int k(0); k < d.Nz; k++) {
          d2.aux[d.id(var, i, j, k)] *= 0.9;
        }
      }
    }
  }
  // Solve and re-check
  model2.getPrimitiveVars(d2.cons, d2.prims, d2.aux);

  for (int var(0); var < d.Nprims; var++) {
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        for (int k(0); k < d.Nz; k++) {
          EXPECT_NEAR(d.prims[d.id(var, i, j, k)], d2.prims[d.id(var, i, j, k)], tol);
        }
      }
    }
  }
  for (int var(0); var < d.Naux; var++) {
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        for (int k(0); k < d.Nz; k++) {
          EXPECT_NEAR(d.aux[d.id(var, i, j, k)], d2.aux[d.id(var, i, j, k)], tol);
        }
      }
    }
  }

}





/* ######################### Test prims2all transform ########################*/

TEST(SRMHD, PrimsToAll)
{
  // Set up
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 1.0, &env);
  SRMHD model(&d);
  Periodic bcs(&d);
  Simulation sim(&d, &env);
  InitialFunc init(&d);

  // Set cons and aux vars to rubbish so we know they have changed, and
  // ser all primitive variables to (nearly) zero: 0 will give zero-division error
  for (int i(0); i < d.Nx; i++) {
    for (int j(0); j < d.Ny; j++) {
      for (int k(0); k < d.Nz; k++) {
        for (int var(0); var < d.Nprims; var++) {
          d.prims[d.id(var, i, j, k)] = 1e-16;
        }
        for (int var(0); var < d.Ncons; var++) {
          d.cons[d.id(var, i, j, k)] = 3.1415926;
        }
        for (int var(0); var < d.Naux; var++) {
          d.aux[d.id(var, i, j, k)] = 2.718281828;
        }
      }
    }
  }

  // Quick check that worked
  EXPECT_EQ(d.cons[d.id(6, 2, 5, 5)], 3.1415926);
  EXPECT_EQ(d.aux[d.id(7, 2, 5, 5)], 2.718281828);
  EXPECT_EQ(d.prims[d.id(7, 8, 4, 5)], 1e-16);

  // Apply conversion and check all cons and aux are zero except h and W
  model.primsToAll(d.cons, d.prims, d.aux);

  for (int i(0); i < d.Nx; i++) {
    for (int j(0); j < d.Ny; j++) {
      for (int k(0); k < d.Nz; k++) {
        // Conserved
        for (int var(0); var < d.Ncons; var++) {
          EXPECT_NEAR(d.cons[d.id(var, i, j, k)], 0.0, 1e-15);
        }
        // Auxiliary
        for (int var(0); var < d.Naux; var++) {
          if (var == 0) EXPECT_EQ(d.aux[d.id(var, i, j, k)], 3.5);
          else if (var == 1) EXPECT_NEAR(d.aux[d.id(var, i, j, k)], 1.0, 1e-15);
          else if (var == 2) EXPECT_EQ(d.aux[d.id(var, i, j, k)], 1.5);
          else if (var == 3) EXPECT_EQ(d.aux[d.id(var, i, j, k)], sqrt(10.0/21.0));
          else EXPECT_NEAR(d.aux[d.id(var, i, j, k)], 0.0, 1e-15);
        }
      }
    }
  }
}
