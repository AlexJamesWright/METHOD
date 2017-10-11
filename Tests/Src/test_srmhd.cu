#include "gtest/gtest.h"
#include "srmhd.h"
#include "simulation.h"
#include "simData.h"
#include "initFunc.h"
#include <cstdlib>
#include <cmath>

namespace {

  /* ##################### Test default model constructor ####################*/

  TEST(SRMHD, DefaultConstructor)
  {
    SRMHD model;
    EXPECT_EQ(model.Ncons, 9);
    EXPECT_EQ(model.Nprims, 8);
    EXPECT_EQ(model.Naux, 10);
  }




  /* ######################### Test model constructor ########################*/

  TEST(SRMHD, Constructor)
  {
    Data d(100, 10, 0, 1, -0.5, 0.5, 0.8);
    SRMHD model(&d);
    EXPECT_EQ(model.data->Ncons, 9);
    EXPECT_EQ(model.data->Nprims, 8);
    EXPECT_EQ(model.data->Naux, 10);

  }



  /* ######################### Test prims2all transform ########################*/

  TEST(SRMHD, PrimsToAll)
  {
    // Set up
    Data d(10, 10, 0, 1, 0, 1, 1.0);
    SRMHD model(&d);
    Simulation sim(&d);
    InitialFunc init(&d);

    // Set cons and aux vars to rubbish so we know they have changed, and
    // ser all primitive variables to (nearly) zero: 0 will give zero-division error
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        for (int var(0); var < d.Nprims; var++) {
          d.prims[d.id(var, i, j)] = 1e-16;
        }
        for (int var(0); var < d.Ncons; var++) {
          d.cons[d.id(var, i, j)] = 3.1415926;
        }
        for (int var(0); var < d.Naux; var++) {
          d.aux[d.id(var, i, j)] = 2.718281828;
        }
      }
    }

    // Quick check that worked
    EXPECT_EQ(d.cons[d.id(6, 2, 5)], 3.1415926);
    EXPECT_EQ(d.aux[d.id(7, 2, 5)], 2.718281828);
    EXPECT_EQ(d.prims[d.id(7, 8, 4)], 1e-16);

    // Apply conversion and check all cons and aux are zero except h and W
    model.primsToAll(d.cons, d.prims, d.aux);

    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        // Conserved
        for (int var(0); var < d.Ncons; var++) {
          EXPECT_NEAR(d.cons[d.id(var, i, j)], 0.0, 1e-15);
        }
        // Auxilliary
        for (int var(0); var < d.Naux; var++) {
          if (var == 0) EXPECT_EQ(d.aux[d.id(var, i, j)], 3.5);
          else if (var == 1) EXPECT_NEAR(d.aux[d.id(var, i, j)], 1.0, 1e-15);
          else if (var == 2) EXPECT_EQ(d.aux[d.id(var, i, j)], 1.5);
          else if (var == 3) EXPECT_EQ(d.aux[d.id(var, i, j)], sqrt(10.0/21.0));
          else EXPECT_NEAR(d.aux[d.id(var, i, j)], 0.0, 1e-15);
        }
      }
    }
  }


  /* ######################### Test source contribution ########################*/

  TEST(SRMHD, SourceTerm)
  {

    // Set up
    Data d(10, 10, 0, 1, 0, 1, 1.0, 0.5, 4, 5.0/3.0, 0.0, 0, 0, 0, 0, 0.5);
    SRMHD model(&d);
    Simulation sim(&d);

    // Set cons to something
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        for (int var(0); var < d.Ncons; var++) {
          d.cons[d.id(var, i, j)] = 3.1415926;
        }
      }
    }
    // Determine source
    model.sourceTerm(d.cons, d.prims, d.aux, d.source);
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        for (int var(0); var < d.Ncons; var++) {
          if (var == 8) EXPECT_EQ(d.source[d.id(var, i, j)], -3.1415926 / (0.25));
          else EXPECT_EQ(d.source[d.id(var, i, j)], 0);
        }
      }
    }
  }


  /* ######################### Test flux vector splitting ########################*/

  TEST(SRMHD, FluxVectorSplittingStationary)
  {

    // Set up
    Data d(10, 10, 0, 1, 0, 1, 1.0, 0.5, 4, 5.0/3.0, 0.0, 0, 0, 0, 0, 0.5);
    SRMHD model(&d);
    Simulation sim(&d);

    // Set state to stationary equilibrium state
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        d.prims[d.id(0, i, j)] = 0.5; // Require non-zero density
        d.prims[d.id(1, i, j)] = 0.0;
        d.prims[d.id(2, i, j)] = 0.0;
        d.prims[d.id(3, i, j)] = 0.0;
        d.prims[d.id(4, i, j)] = 0.0;
        d.prims[d.id(5, i, j)] = 0.0;
        d.prims[d.id(6, i, j)] = 0.0;
        d.prims[d.id(7, i, j)] = 0.0;
        d.prims[d.id(8, i, j)] = 0.0;
      }
    }

    model.primsToAll(d.cons, d.prims, d.aux);

    // System is stationary, there should be zero flux
    // x-direction
    model.fluxFunc(d.cons, d.prims, d.aux, d.f, d.fnet, 0);
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        for (int var(0); var < d.Ncons; var++) {
          EXPECT_EQ(d.fnet[d.id(var, i, j)], 0.0);
        }
      }
    }
    // y-direction
    model.fluxFunc(d.cons, d.prims, d.aux, d.f, d.fnet, 1);
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        for (int var(0); var < d.Ncons; var++) {
          EXPECT_EQ(d.fnet[d.id(var, i, j)], 0.0);
        }
      }
    }
  }


}
